#include "HillClimbers_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

#include <thrust/sequence.h>

#define MAX_MORTALITY 0.975

void fecundity_genotype_phenotype_map::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = raw_pointer_cast(&species->fgenotype[0][0]);
	float *fgen2 = raw_pointer_cast(&species->fgenotype[1][0]);
	float *fgen3 = raw_pointer_cast(&species->fgenotype[2][0]);
	float *fgen4 = raw_pointer_cast(&species->fgenotype[3][0]);
	float *fgen5 = raw_pointer_cast(&species->fgenotype[4][0]);
	
	float *mgen1 = raw_pointer_cast(&species->mgenotype[0][0]);
	float *mgen2 = raw_pointer_cast(&species->mgenotype[1][0]);
	float *mgen3 = raw_pointer_cast(&species->mgenotype[2][0]);
	float *mgen4 = raw_pointer_cast(&species->mgenotype[3][0]);
	float *mgen5 = raw_pointer_cast(&species->mgenotype[4][0]);
	
	// wrap the allelic_effects coefficient parameters from vectors to arrays:
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
	float *coefficient_1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *coefficient_2 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *coefficient_3 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *coefficient_4 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF4")[0]);

	//Set up functor
	fecundity_calculator fecundity_calculator_functor(fgen1, fgen2, fgen3, fgen4, fgen5, mgen1, mgen2, mgen3, mgen4, mgen5, coefficient_0, coefficient_1, coefficient_2, coefficient_3, coefficient_4);

	thrust::device_vector<int> individuals_indices(num_kids);
	thrust::sequence(individuals_indices.begin(), individuals_indices.begin() + num_kids, index_case);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
					   individuals_indices.begin(), 
					   species->phenotype[phenotype_index].begin() + index_case,
					   species->deme.begin() + index_case)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
					   individuals_indices.begin() + num_kids, 
					   species->phenotype[phenotype_index].begin() + index_case + num_kids,
					   species->deme.begin() + index_case + num_kids)),
			fecundity_calculator_functor);
	}

struct max_mortality
  {
    __host__ __device__
    bool operator()(float x)
    {
      return x > MAX_MORTALITY;
    }
  };

void mortality_genotype_phenotype_map::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = raw_pointer_cast(&species->fgenotype[5][0]);
	float *fgen2 = raw_pointer_cast(&species->fgenotype[6][0]);
	float *fgen3 = raw_pointer_cast(&species->fgenotype[7][0]);
	float *fgen4 = raw_pointer_cast(&species->fgenotype[8][0]);
	float *fgen5 = raw_pointer_cast(&species->fgenotype[9][0]);
	
	float *mgen1 = raw_pointer_cast(&species->mgenotype[5][0]);
	float *mgen2 = raw_pointer_cast(&species->mgenotype[6][0]);
	float *mgen3 = raw_pointer_cast(&species->mgenotype[7][0]);
	float *mgen4 = raw_pointer_cast(&species->mgenotype[8][0]);
	float *mgen5 = raw_pointer_cast(&species->mgenotype[9][0]);
	
	// wrap the allelic_effects coefficient parameters from vectors to arrays:
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
	float *coefficient_1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *coefficient_2 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *coefficient_3 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *coefficient_4 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF4")[0]);

	//Set up functor
	mortality_calculator mortality_calculator_functor(fgen1, fgen2, fgen3, fgen4, fgen5, mgen1, mgen2, mgen3, mgen4, mgen5, coefficient_0, coefficient_1, coefficient_2, coefficient_3, coefficient_4);

	thrust::device_vector<int> individuals_indices(num_kids);
	thrust::sequence(individuals_indices.begin(), individuals_indices.begin() + num_kids, index_case);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
					   individuals_indices.begin(), 
					   species->phenotype[phenotype_index].begin() + index_case,
					   species->deme.begin() + index_case)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
					   individuals_indices.begin() + num_kids, 
					   species->phenotype[phenotype_index].begin() + index_case + num_kids,
					   species->deme.begin() + index_case + num_kids)),
			mortality_calculator_functor);

	max_mortality predicate;
	
	thrust::replace_if(species->phenotype[phenotype_index].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case, predicate, MAX_MORTALITY);
	}
 
