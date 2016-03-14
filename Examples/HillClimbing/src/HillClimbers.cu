#include "HillClimbers.h"

#include <species/add_kids/genotype_phenotype_map.h>
#include <species/update/updatebehavior.h>

#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>

HillClimbers::HillClimbers(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds_stochastic(size_val, maxsize_val, seed_val, ndemes, species_ID_val)
	{
	// Assume everyone starts at age = 0.
	thrust::fill(age.begin(), age.begin() + size, 0);

	initialize_demes();
	
	// Specify the genetics by assuming allelic values are gaussian-distributed
	for (int i=0; i < nloci; i++)
		{		
		draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
		draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
		}

	//Set phenotype
	setPhenotype(0, size);

	assignSex(0, size);
	}

void HillClimbers::addKids()
	{   
	demeCalculations();

	HillClimberParents *exampleParents;
	exampleParents = new HillClimberParents(this);
	exampleParents->setup_parents();

	if (exampleParents->Potential_Number_of_Kids > 0)
		{
		EggsNeonates *	HillClimber_Neonates;
		HillClimber_Neonates = new EggsNeonates (this, exampleParents->kids_per_mom);
		HillClimber_Neonates->inherit_genotypes(exampleParents->probability_individual_becomes_female_parent,  exampleParents->probability_individual_becomes_male_parent);

		setPhenotype(HillClimber_Neonates->previous_pop_size, HillClimber_Neonates->Total_Number_of_Neonates);
		assignSex(HillClimber_Neonates->previous_pop_size, HillClimber_Neonates->Total_Number_of_Neonates);
		delete HillClimber_Neonates;     
		sortByDeme();
		}
	demeCalculations();
	delete exampleParents;
	}

void HillClimbers::setPhenotype(int index, int num_inds_to_calculate)
	{   
	for (int i=0; i < nphen; i++)
		{   
		GenotypePhenotypeMap *genotype_phenotype_map;
		genotype_phenotype_map = genotype_phenotype_map->create_genotype_phenotype_map(this, i, index, num_inds_to_calculate);
		genotype_phenotype_map->calculate_phenotype( this );
		delete genotype_phenotype_map;
		}
	}

void HillClimbers::assignSex(int index, int num_inds_to_calculate)
	{
	thrust::host_vector<int> neonates_sex(num_inds_to_calculate);
	draw_bernoulli(num_inds_to_calculate, 0.5, neonates_sex, gen);
	thrust::copy(neonates_sex.begin(), neonates_sex.begin() + num_inds_to_calculate, sex.begin() + index);
	}

void HillClimbers::update(inds_stochastic **species)
	{
	UpdateBehavior *updatebehavior;
	updatebehavior = updatebehavior->create_updateBehavior(species, NULL, 0);
	updatebehavior->update();
	}


void HillClimbers::initialize_demes()
	{
	thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);

	for (int i=0; i < Num_Demes; i++)
		{
		float startsize = (float) (size/Num_Demes);
		int temp1 = (int) startsize * i;
		int temp2 = (int) startsize * (i+1);
		thrust::fill(deme.begin() + temp1, deme.begin() + temp2, i);
		}
	}


