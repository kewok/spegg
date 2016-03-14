#include <thrust/transform.h>
#include <thrust/functional.h>

#include "coevolvingSpeciesParents.h"
#include "coevolvingSpecie.h"
#include "coevolvingSpecies_mating_kernel_functors.h"

coevolvingSpeciesParents::coevolvingSpeciesParents(coevolvingSpecie *species) : Parents(species)
	{
	FECUNDITY_PHENOTYPE_INDEX = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	age.resize(size);
	thrust::copy(species->age.begin(), species->age.begin() + size, age.begin());
	}

void coevolvingSpeciesParents::determine_male_parent_eligibility()
	{
	thrust::fill(will_reproduceM.begin(), will_reproduceM.begin() + size, 0);
	
	float *m_ages_at_maturity_ptr = demeParameters->get_vector_ptr("Age_at_maturity");

	if (demeParameters->species_specific_values["reproduction_mode"]==0.0)
		{	
		bisexual_parental_eligibility_functor bisexual_parental_eligibility(m_ages_at_maturity_ptr);

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceM.begin(), 
							   deme.begin(), 
							   age.begin())),
		 		 thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceM.begin() + size, 
							   deme.begin() + size, 
							   age.begin() + size)),
			 	 bisexual_parental_eligibility);
		}
	if (demeParameters->species_specific_values["reproduction_mode"]==1.0)
		{
		coevolving_parental_eligibility_functor parental_eligibility(m_ages_at_maturity_ptr, 1);

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceM.begin(),
							   sex.begin(), 
							   deme.begin(), 
							   age.begin())),
		 		 thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceM.begin() + size, 
							   sex.begin(),
							   deme.begin() + size, 
							   age.begin() + size)),
			 	 parental_eligibility);
		}
	}

void coevolvingSpeciesParents::determine_female_parent_eligibility()
	{
	thrust::fill(will_reproduceF.begin(), will_reproduceF.begin() + size, 0);
	
	float *f_ages_at_maturity_ptr = demeParameters->get_vector_ptr("Age_at_maturity");

	if (demeParameters->species_specific_values["reproduction_mode"]==0.0)
		{	
		bisexual_parental_eligibility_functor bisexual_parental_eligibility(f_ages_at_maturity_ptr);

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceF.begin(), 
							   deme.begin(), 
							   age.begin())),
		 		 thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceF.begin() + size, 
							   deme.begin() + size, 
							   age.begin() + size)),
			 	 bisexual_parental_eligibility);
		}
	if (demeParameters->species_specific_values["reproduction_mode"]==1.0)
		{
		coevolving_parental_eligibility_functor parental_eligibility(f_ages_at_maturity_ptr, 0);

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceF.begin(),
							   sex.begin(), 
							   deme.begin(), 
							   age.begin())),
		 		 thrust::make_zip_iterator(thrust::make_tuple(
							   will_reproduceF.begin() + size, 
							   sex.begin(),
							   deme.begin() + size, 
							   age.begin() + size)),
			 	 parental_eligibility);
		}
	}

void coevolvingSpeciesParents::determine_probability_individual_becomes_female_parent()
	{
	thrust::multiplies<float> op;
	thrust::transform(will_reproduceF.begin(), will_reproduceF.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),  probability_individual_becomes_female_parent.begin(), op);
	}

void coevolvingSpeciesParents::determine_probability_individual_becomes_male_parent()
	{
	thrust::multiplies<float> op;
	thrust::transform(will_reproduceM.begin(), will_reproduceM.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),  probability_individual_becomes_male_parent.begin(), op);
	}

