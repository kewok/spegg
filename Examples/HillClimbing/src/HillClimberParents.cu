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
	int deme_offset = 0;
	float female_reproductive_skew = *(demeParameters->get_vector_ptr("FEMALE_REPRODUCTIVE_SKEW") + deme_offset);

	for (int i=0; i < size; i++)
		{
		// added to minimize memory access:
		if (deme[i] != deme_offset)
			{
			deme_offset = deme[i];
			female_reproductive_skew = *(demeParameters->get_vector_ptr("FEMALE_REPRODUCTIVE_SKEW") + deme_offset);
			}
		if (will_reproduceF[i] > 0)
			{
			probability_individual_becomes_female_parent[i] = powf(phenotype[FECUNDITY_PHENOTYPE_INDEX][i], female_reproductive_skew);
			}
		}
	}

void HillClimberParents::determine_probability_individual_becomes_male_parent()
	{
	int deme_offset = 0;
	float male_reproductive_skew = *(demeParameters->get_vector_ptr("MALE_REPRODUCTIVE_SKEW") + deme_offset);

	for (int i=0; i < size; i++)
		{
		// added to minimize memory access:
		if (deme[i] != deme_offset)
			{
			deme_offset = deme[i];
			male_reproductive_skew = *(demeParameters->get_vector_ptr("MALE_REPRODUCTIVE_SKEW") + deme_offset);
			}
		if (will_reproduceM[i] > 0)
			{
			probability_individual_becomes_male_parent[i] = powf(phenotype[FECUNDITY_PHENOTYPE_INDEX][i], male_reproductive_skew);
			}
		}
	}

