#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>

#include "HillClimberParents.h"
#include "HillClimbers.h"

HillClimberParents::HillClimberParents(HillClimbers *species) : Parents(species)
	{
	FECUNDITY_PHENOTYPE_INDEX = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	}
	
void HillClimberParents::determine_probability_individual_becomes_female_parent()
	{
	float *reproductive_skews = raw_pointer_cast(&demeParameters->get_vector_ptr("FEMALE_REPRODUCTIVE_SKEW")[0]);

	reproductive_probability_functor reproductive_probability_calculator(reproductive_skews);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), deme.begin(), phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(), probability_individual_becomes_female_parent.begin())),
		 	 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, deme.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + size, probability_individual_becomes_female_parent.begin() + size)),
			 reproductive_probability_calculator);	
	}

void HillClimberParents::determine_probability_individual_becomes_male_parent()
	{
	float *reproductive_skews = raw_pointer_cast(&demeParameters->get_vector_ptr("MALE_REPRODUCTIVE_SKEW")[0]);

	reproductive_probability_functor reproductive_probability_calculator(reproductive_skews);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin(), deme.begin(), phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(), probability_individual_becomes_male_parent.begin())),
		 	 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin() + size, deme.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + size, probability_individual_becomes_male_parent.begin() + size)),
			 reproductive_probability_calculator);	
	}

