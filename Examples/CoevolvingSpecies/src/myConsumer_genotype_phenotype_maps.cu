#include <thrust/gather.h>

#include "myConsumer_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

void consumer_fecundity_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void consumer_attack_phenotype::calculate_phenotype(inds *species)
	{
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
	float *coefficient_1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *coefficient_2 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);

	consumer_attack_phenotype_calculator consumer_attack_functor(coefficient_0, coefficient_1, coefficient_2);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
						 species->deme.begin() + index_case, 
						 species->fgenotype[0].begin() + index_case, 
						 species->mgenotype[0].begin() + index_case, 
						 species->fgenotype[1].begin() + index_case, 
						 species->mgenotype[1].begin() + index_case, 
						 species->fgenotype[2].begin() + index_case, 
						 species->mgenotype[2].begin() + index_case, 
						 species->phenotype[CONSUMER_ATTACK_PHENOTYPE_INDEX].begin() + index_case)),
			thrust::make_zip_iterator(thrust::make_tuple(
						  species->deme.begin() + index_case + num_kids,
						  species->fgenotype[0].begin() + index_case + num_kids,
						  species->mgenotype[0].begin() + index_case + num_kids,
						  species->fgenotype[1].begin() + index_case + num_kids,
						  species->mgenotype[1].begin() + index_case + num_kids,
						  species->fgenotype[2].begin() + index_case + num_kids,
						  species->mgenotype[2].begin() + index_case + num_kids,
						  species->phenotype[CONSUMER_ATTACK_PHENOTYPE_INDEX].begin() + index_case + num_kids)),
			consumer_attack_functor);
	}


void consumer_mortality_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[MORTALITY_PHENOTYPE_INDEX].begin() + index_case);
	}
