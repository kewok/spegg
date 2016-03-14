#ifndef COEVOLVING_SPECIES_PARENTS_CLASS_H
#define COEVOLVING_SPECIES_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>

class coevolvingSpeciesParents : public Parents
	{
	class coevolvingSpecie *species;
	public:
		coevolvingSpeciesParents(coevolvingSpecie *species);

	protected:
		int FECUNDITY_PHENOTYPE_INDEX;
		thrust::device_vector<int> age;

		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();

		void determine_female_parent_eligibility();
		void determine_male_parent_eligibility();
	};


#endif
