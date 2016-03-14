#ifndef HILLCLIMBER_PARENTS_CLASS_H
#define HILLCLIMBER_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>

class HillClimberParents : public Parents
	{
	class HillClimbers *species;
	public:
		HillClimberParents(HillClimbers *species);
		
	protected:
		int FECUNDITY_PHENOTYPE_INDEX;

		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();
	};

#endif 
