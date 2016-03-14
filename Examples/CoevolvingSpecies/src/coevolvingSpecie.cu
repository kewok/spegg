#include "coevolvingSpecie.h"

#include <species/add_kids/genotype_phenotype_map.h>
#include <species/update/updatebehavior.h>
#include <species/movement/migration_behavior.h>

#include <util/footimer2.h>

#include "myResource.h"
#include "myConsumer.h"

#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>

coevolvingSpecie::coevolvingSpecie(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds_stochastic_migratory(size_val, maxsize_val, seed_val, ndemes, species_ID_val)
	{
	species_ID = species_ID_val;
	// Assume everyone starts at age = 0.
	thrust::fill(age.begin(), age.begin() + size, 0);
	demeParameters = new DemeSettings_coevolvingSpecies("deme_config.txt", species_ID);
	}

void coevolvingSpecie::addKids()	
	{
	demeCalculations();

	coevolvingSpeciesParents *exampleParents;
	exampleParents = new coevolvingSpeciesParents(this);
	exampleParents->setup_parents();

	if (exampleParents->Potential_Number_of_Kids > 0)
		{
		EggsNeonates *exampleNeonates;
		exampleNeonates = new EggsNeonates (this, exampleParents->kids_per_mom);
		exampleNeonates->inherit_genotypes(exampleParents->probability_individual_becomes_female_parent,  exampleParents->probability_individual_becomes_male_parent);
		setPhenotype(exampleNeonates->previous_pop_size, exampleNeonates->Total_Number_of_Neonates);
		assignSex(exampleNeonates->previous_pop_size, exampleNeonates->Total_Number_of_Neonates);
		delete exampleNeonates;     
		sortByDeme();
		}

	demeCalculations();
	delete exampleParents;
	}


void coevolvingSpecie::update(inds_stochastic **species)
	{
	UpdateBehavior *updatebehavior;
	updatebehavior = updatebehavior->create_updateBehavior(species, NULL, species_ID);
	updatebehavior->update();
	updatebehavior->determine_mortality(this);
	delete updatebehavior;

	// Increment ages by one.
	thrust::host_vector<int> ones(size);
	thrust::fill(ones.begin(), ones.end(), 1);
	thrust::transform_if(age.begin(), age.begin() + size, status.begin(), age.begin(), unary_plus<int>(1), thrust::identity<int>());

	// Update the migration propensities; note that if we model the individual migration rates as an evolvable phenotype, then this code might look something like:
	/*
	thrust::copy(phenotype[MIGRATION_PHENOTYPE_INDEX]->begin(), phenotype[MIGRATION_PHENOTYPE_INDEX]->begin() + size, individual_migration_rate.begin());
	*/
	individual_migration_rate.resize(size);
	thrust::gather(deme.begin(), deme.begin() + size, demeParameters->get_vector_ptr("individual_migration_propensity"), individual_migration_rate.begin());
	}

void coevolvingSpecie::setPhenotype(int index, int num_inds_to_calculate)
	{	
	for (int i=0; i < nphen; i++)
		{
		GenotypePhenotypeMap *genotype_phenotype_map;
		genotype_phenotype_map = genotype_phenotype_map->create_genotype_phenotype_map(this, i, index, num_inds_to_calculate);
		genotype_phenotype_map->calculate_phenotype(this);
		delete genotype_phenotype_map;
		}
	}

void coevolvingSpecie::initialize_demes()
	{
	// divide up the demes; for now, just assume all patches share equal starting sizes (sizes/Num_Demes)
	thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);

	for (int i=0; i < Num_Demes; i++)
		{
		float startsize = (float) (size/Num_Demes);
		int temp1 = (int) startsize * i;
		int temp2 = (int) startsize * (i+1);
		thrust::fill(deme.begin() + temp1, deme.begin() + temp2, i);
		}
	// Set the deme sizes accordingly
	demeCalculations();
	}

void coevolvingSpecie::assignSex(int index, int n)
	{
	thrust::host_vector<int> neonates_sex(n);
	draw_bernoulli(n, 0.5, neonates_sex, gen);
	thrust::copy(neonates_sex.begin(), neonates_sex.begin() + n, sex.begin() + index);
	}

void coevolvingSpecie::migrate()
	{
	MigrationBehavior *migration_behavior;
	migration_behavior = migration_behavior->create_migrationBehavior(this);
	migration_behavior->migrate();
	delete migration_behavior;
	}


