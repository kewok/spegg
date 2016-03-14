#include <thrust/gather.h>

#include "myConsumer_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

void consumer_fecundity_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void consumer_attack_phenotype::calculate_phenotype(inds *species)
	{
	float *coefficient_0 = &Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0];
	float *coefficient_1 = &Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0];
	float *coefficient_2 = &Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0];

	int deme_offset = 0;
	float COEF_0 = *(coefficient_0 + deme_offset);
	float COEF_1 = *(coefficient_1 + deme_offset);
	float COEF_2 = *(coefficient_2 + deme_offset);

	for (int i=index_case; i < index_case + num_kids; i++)
		{
		if (species->deme[i] != deme_offset)
			{
			deme_offset = species->deme[i];

			COEF_0 = *(coefficient_0 + deme_offset);
			COEF_1 = *(coefficient_1 + deme_offset);
			COEF_2 = *(coefficient_2 + deme_offset);
			}
		species->phenotype[CONSUMER_ATTACK_PHENOTYPE_INDEX][i] = COEF_0*(species->fgenotype[0][i] + species->mgenotype[0][i]) + COEF_1*(species->fgenotype[1][i] + species->mgenotype[1][i]) + COEF_2*(species->fgenotype[2][i] + species->mgenotype[2][i]);	
		species->phenotype[CONSUMER_ATTACK_PHENOTYPE_INDEX][i] = species->phenotype[CONSUMER_ATTACK_PHENOTYPE_INDEX][i]*0.5;
		}
	}


void consumer_mortality_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[MORTALITY_PHENOTYPE_INDEX].begin() + index_case);
	}
