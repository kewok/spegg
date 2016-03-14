#ifndef PENGUIN_PARENTS_CLASS_H
#define PENGUIN_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>

class Penguin_Parents : public Parents
	{
	class Penguins *species;
	public:
		Penguin_Parents(Penguins *species);
		
	protected:
		int FECUNDITY_PHENOTYPE_INDEX;
		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();
	};

#endif
