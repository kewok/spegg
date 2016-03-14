#ifndef FISH_PARENTS_CLASS_H
#define FISH_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>
#include "Fish_mating_kernel_functors.h"

class Fish_Parents : public Parents
	{
	class Fish *species;
	public:
		Fish_Parents(Fish *species);
		
	protected:
		int IRREVERSIBLE_MASS_PHENOTYPE;
		int FECUNDITY_PHENOTYPE;
		int REVERSIBLE_MASS_PHENOTYPE;

		void determine_female_parent_eligibility();
		void determine_male_parent_eligibility();
		void female_fecundity();

		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();
		void female_nesting_success();		
	};

#endif 
