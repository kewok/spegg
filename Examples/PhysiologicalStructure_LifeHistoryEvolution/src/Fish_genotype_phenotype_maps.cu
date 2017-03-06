#include "Fish_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>

void fecundity_phenotype::calculate_phenotype(inds *species)
	{
	thrust::fill(species->phenotype[phenotype_index].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case + num_kids, 0.0);
	}


void mortality_phenotype::calculate_phenotype(inds *species)
	{	
	thrust::fill(species->phenotype[phenotype_index].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case + num_kids, 0.0);
	}

void irreversible_mass_at_birth::calculate_phenotype(inds *species)
	{
	int irreversible_mass_phenotype_index = species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];

	// For multiple loci, use:	for (int i=0; i < Parameters-> get_Number_of_Loci(); i++)
	// wrap the genotypes:
	float *fgen1 = (&species->fgenotype[1][0]);
	float *fgen2 = (&species->fgenotype[2][0]);
	float *fgen3 = (&species->fgenotype[3][0]);
	float *fgen4 = (&species->fgenotype[4][0]);
	float *fgen5 = (&species->fgenotype[5][0]);
	float *fgen6 = (&species->fgenotype[6][0]);
	float *fgen7 = (&species->fgenotype[7][0]);

	float *mgen1 = (&species->mgenotype[1][0]);
	float *mgen2 = (&species->mgenotype[2][0]);
	float *mgen3 = (&species->mgenotype[3][0]);
	float *mgen4 = (&species->mgenotype[4][0]);
	float *mgen5 = (&species->mgenotype[5][0]);
	float *mgen6 = (&species->mgenotype[6][0]);
	float *mgen7 = (&species->mgenotype[7][0]);

// wrap the allelic effects:
	float *locus_effect1 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF0")[0]);
	float *locus_effect2 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF1")[0]);
	float *locus_effect3 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF2")[0]);
	float *locus_effect4 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF3")[0]);
	float *locus_effect5 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF4")[0]);
	float *locus_effect6 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF5")[0]);
	float *locus_effect7 = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF6")[0]);
	float *epistatic_effect = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_GENPHEN_MAP_COEF7")[0]);
	float *intercept_ptr = (&Parameters->get_vector_ptr("OFFSPRING_SIZE_CONSTANT")[0]);
	float *maternal_effects_ptr = (&species->phenotype[irreversible_mass_phenotype_index][0]);

	offspring_size_calculator offspring_size_functor(intercept_ptr, fgen1, fgen2, fgen3, fgen4, fgen5, fgen6, fgen7, mgen1, mgen2, mgen3, mgen4, mgen5, mgen6, mgen7, locus_effect1, locus_effect2, locus_effect3, locus_effect4, locus_effect5, locus_effect6, locus_effect7, epistatic_effect, maternal_effects_ptr);

	// draw the gaussian for phenotype deviation and create list of individuals
	thrust::host_vector<float> deviation(num_kids);

	for (int i=0; i < num_kids; i++)
		{
		deviation[i] = gsl_ran_gaussian(gen, 0.0);
		}
	
	thrust::host_vector<int> individuals(index_case + num_kids);
	thrust::sequence(individuals.begin(), individuals.begin() + index_case + num_kids, 0);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + index_case, deviation.begin(), species->phenotype[phenotype_index].begin() + index_case, species->deme.begin() + index_case, species->maternal_id.begin() + index_case)),
	                 thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + index_case + num_kids, deviation.begin() + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids, species->deme.begin() + index_case + num_kids, species->maternal_id.begin() + index_case + num_kids)),
	                 offspring_size_functor);

	thrust::copy(species->phenotype[phenotype_index].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case + num_kids, species->phenotype[irreversible_mass_phenotype_index].begin() + index_case);
	}

void reversible_mass_at_birth::calculate_phenotype(inds *species)
	{
	int irreversible_mass_phenotype_index = species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];

	thrust::host_vector<float> starting_conditions(num_kids);
	float *juvenile_maximum_condition_ptr = species->demeParameters->get_vector_ptr("juvenile_maximum_condition");
	for (int i=0; i < num_kids; i++)
		{
		// Make sure this is working as intended:
		starting_conditions[i] = *(juvenile_maximum_condition_ptr + species->deme[index_case + i]);
		species->phenotype[phenotype_index][i + index_case] = starting_conditions[i] * species->phenotype[irreversible_mass_phenotype_index][i + index_case];
		}
	}

void satiation_at_birth::calculate_phenotype(inds *species )
	{
	thrust::fill(species->phenotype[phenotype_index].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case + num_kids, 0);
	}

void eggsize_at_birth::calculate_phenotype(inds *species )
	{
	int irreversible_mass_phenotype_index = species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];

	thrust::copy(species->phenotype[irreversible_mass_phenotype_index].begin() + index_case, species->phenotype[irreversible_mass_phenotype_index].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case);
	}

