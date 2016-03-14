#include "Penguins_genotype_phenotype_maps.h"

void fecundity_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
	float *coefficient_1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);

	//Instantiate functor
	fecundity_calculator fecundity_calculator_functor(constants, coefficient_0, coefficient_1);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->fgenotype[1].begin() + index_case, species->mgenotype[1].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)),
			thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->fgenotype[1].begin() + index_case + num_kids, species->mgenotype[1].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)),
			fecundity_calculator_functor);
	}


void mortality_phenotype::calculate_phenotype(inds *species )
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

	//Set up functor
	mortality_calculator mortality_calculator_functor(constants, coefficient_0);

	//Perform mortality operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)),
	                 thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)),
	                 mortality_calculator_functor);
	}


void crown_color_phenotype::calculate_phenotype(inds *species )
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

	//Set up functor
	crown_color_calculator crown_color_functor(constants, coefficient_0);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[2].begin() + index_case, species->mgenotype[2].begin() + index_case, species->phenotype[CROWN_COLOR_INDEX].begin() + index_case)),
	                 thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[2].begin() + index_case + num_kids, species->mgenotype[2].begin() + index_case + num_kids, species->phenotype[CROWN_COLOR_INDEX].begin() + index_case + num_kids)),
	                 crown_color_functor);
	}
