#include "Fish_Parents.h"
#include "Fish.h"
#include "Fish_mating_kernel_functors.h"

Fish_Parents::Fish_Parents(Fish *species) : Parents(species)
	{
	this->IRREVERSIBLE_MASS_PHENOTYPE = demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];
	this->FECUNDITY_PHENOTYPE = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE"];
	this->REVERSIBLE_MASS_PHENOTYPE = demeParameters->species_specific_values["REVERSIBLE_MASS_PHENOTYPE"];
	}
	
void Fish_Parents::determine_female_parent_eligibility()
	{
	thrust::fill(will_reproduceF.begin(), will_reproduceF.begin() + size, 0);
	
	//Set up eligibility functor.
	float *f_sizes_at_maturity_ptr = raw_pointer_cast(demeParameters->get_vector_ptr("F_sizes_at_maturity"));
	parental_eligibility_functor parental_eligibility(f_sizes_at_maturity_ptr, 0);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), sex.begin(), deme.begin(), phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin())),
		 	thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, sex.begin() + size, deme.begin() + size, phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size)),
			 parental_eligibility);	
	}

void Fish_Parents::determine_male_parent_eligibility()
	{
	thrust::fill(will_reproduceM.begin(), will_reproduceM.begin() + size, 0);
	
	//Set up eligibility functor.
	float *m_sizes_at_maturity_ptr = raw_pointer_cast(demeParameters->get_vector_ptr("M_sizes_at_maturity"));
	parental_eligibility_functor parental_eligibility(m_sizes_at_maturity_ptr, 1);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin(), sex.begin(), deme.begin(), phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin() + size, sex.begin() + size, deme.begin() + size, phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size)),
			 parental_eligibility);	
	}

void Fish_Parents::female_fecundity()
	{
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), phenotype[FECUNDITY_PHENOTYPE].begin(), probability_individual_becomes_female_parent.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, phenotype[FECUNDITY_PHENOTYPE].begin() + size, probability_individual_becomes_female_parent.begin() + size)),
			 female_fecundity_functor());

	thrust::copy(probability_individual_becomes_female_parent.begin(), probability_individual_becomes_female_parent.end(), kids_per_mom.begin());
	}

void Fish_Parents::determine_probability_individual_becomes_female_parent()
	{
	float *f_reproductive_advantage = raw_pointer_cast(demeParameters->get_vector_ptr("F_reproductive_advantage"));
	reproductive_probability_functor F_reproductive_probability(f_reproductive_advantage, 0);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), deme.begin(), phenotype[FECUNDITY_PHENOTYPE].begin(), probability_individual_becomes_female_parent.begin())),
        	 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, deme.begin() + size, phenotype[FECUNDITY_PHENOTYPE].begin() + size, probability_individual_becomes_female_parent.begin() + size)),
               	 F_reproductive_probability);
	}

void Fish_Parents::determine_probability_individual_becomes_male_parent()
	{
	//Set up reproductive probability functor.
	float *m_reproductive_advantage = raw_pointer_cast(demeParameters->get_vector_ptr("M_reproductive_advantage"));

	reproductive_probability_functor M_reproductive_probability(m_reproductive_advantage, 1);

	// For now, assume male reproductive success just depends on his irreversible mass which is proportional to body length 

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin(), deme.begin(), phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin(), probability_individual_becomes_male_parent.begin())),
        	 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin() + size, deme.begin() + size,phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size, probability_individual_becomes_male_parent.begin() + size)),
               	 M_reproductive_probability);
	}

