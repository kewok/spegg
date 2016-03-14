#include <thrust/gather.h>
#include <thrust/sequence.h>

#include "myResource_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

void resource_fecundity_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void resource_mortality_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[MORTALITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void resource_defense_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = raw_pointer_cast(&species->fgenotype[0][0]);
	float *fgen2 = raw_pointer_cast(&species->fgenotype[1][0]);
	float *fgen3 = raw_pointer_cast(&species->fgenotype[2][0]);
	float *fgen9 = raw_pointer_cast(&species->fgenotype[8][0]);
	float *fgen10 = raw_pointer_cast(&species->fgenotype[9][0]);

	float *mgen1 = raw_pointer_cast(&species->mgenotype[0][0]);
	float *mgen2 = raw_pointer_cast(&species->mgenotype[1][0]);
	float *mgen3 = raw_pointer_cast(&species->mgenotype[2][0]);
	float *mgen9 = raw_pointer_cast(&species->mgenotype[8][0]);
	float *mgen10 = raw_pointer_cast(&species->mgenotype[9][0]);

// wrap the allelic effects:
	float *locus_effect1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *locus_effect2 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *locus_effect3 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *locus_effect9 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF9")[0]);
	float *locus_effect10 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF10")[0]);

	resource_defense_calculator resource_defense_functor(fgen1, fgen2, fgen3, fgen9, fgen10, mgen1, mgen2, mgen3, mgen9, mgen10, locus_effect1, locus_effect2, locus_effect3, locus_effect9, locus_effect10);

	thrust::device_vector<int> individuals(index_case + num_kids);
	thrust::sequence(individuals.begin(), individuals.begin() + index_case + num_kids, 0);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
						   individuals.begin() + index_case,
						   species->phenotype[RESOURCE_DEFENSE_PHENOTYPE_INDEX].begin() + index_case,
						   species->deme.begin() + index_case)),		
			 thrust::make_zip_iterator(thrust::make_tuple(
						   individuals.begin() + index_case + num_kids,
						   species->phenotype[RESOURCE_DEFENSE_PHENOTYPE_INDEX].begin() + index_case + num_kids,
						   species->deme.begin() + index_case  + num_kids)),
			 resource_defense_functor);
	}

void resource_competition_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = raw_pointer_cast(&species->fgenotype[0][0]);
	float *fgen2 = raw_pointer_cast(&species->fgenotype[1][0]);
	float *fgen3 = raw_pointer_cast(&species->fgenotype[2][0]);
	float *fgen8 = raw_pointer_cast(&species->fgenotype[7][0]);
	float *fgen10 = raw_pointer_cast(&species->fgenotype[9][0]);

	float *mgen1 = raw_pointer_cast(&species->mgenotype[0][0]);
	float *mgen2 = raw_pointer_cast(&species->mgenotype[1][0]);
	float *mgen3 = raw_pointer_cast(&species->mgenotype[2][0]);
	float *mgen8 = raw_pointer_cast(&species->mgenotype[7][0]);
	float *mgen10 = raw_pointer_cast(&species->mgenotype[9][0]);

// wrap the allelic effects:
	float *locus_effect1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *locus_effect2 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *locus_effect3 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *locus_effect8 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF8")[0]);
	float *locus_effect10 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF10")[0]);

	resource_competitive_ability_calculator resource_competitive_ability_functor(fgen1, fgen2, fgen3, fgen8, fgen10, mgen1, mgen2, mgen3, mgen8, mgen10, locus_effect1, locus_effect2, locus_effect3, locus_effect8, locus_effect10);

	thrust::device_vector<int> individuals(index_case + num_kids);
	thrust::sequence(individuals.begin(), individuals.begin() + index_case + num_kids, 0);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
						   individuals.begin() + index_case,
						   species->phenotype[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX].begin() + index_case,
						   species->deme.begin() + index_case)),		
			 thrust::make_zip_iterator(thrust::make_tuple(
						   individuals.begin() + index_case + num_kids,
						   species->phenotype[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX].begin() + index_case + num_kids,
						   species->deme.begin() + index_case  + num_kids)),
			 resource_competitive_ability_functor);
	}
