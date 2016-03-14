#ifndef UPDATE_COEVOLVING_SPECIES_H
#define UPDATE_COEVOLVING_SPECIES_H

#include <species/update/updatebehavior.h>
#include "coevolvingSpecie.h"
#include "interaction_functor.h"

class update_coevolvingSpecies : public UpdateBehavior
	{
	class coevolvingSpecie **species;
	// Constructor
	public:
		update_coevolvingSpecies(inds_stochastic **species, int species_ID) 
		 	{
			//coevolvingSpecie *specie_ptr = reinterpret_cast <coevolvingSpecie*> (*species);
			this->species = (coevolvingSpecie **) species;
			this->FOCAL_SPECIES_INDEX = species_ID;
		
			// Copy the constants 
			this->size = species[FOCAL_SPECIES_INDEX]->size;
			this->Number_of_Demes = species[FOCAL_SPECIES_INDEX]->Num_Demes;


			this->FECUNDITY_PHENOTYPE_INDEX = species[FOCAL_SPECIES_INDEX]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
			this->MORTALITY_PHENOTYPE_INDEX = species[FOCAL_SPECIES_INDEX]->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];	
			}

		void update();

	protected:
		int FOCAL_SPECIES_INDEX;
		int size;
		int Number_of_Demes;
		int ALTERNATIVE_SPECIES_INDEX;
		int number_of_species;

		thrust::host_vector<int> cumulative_deme_sizes;
		int *cumulative_deme_sizes_ptr;
		
		int phenotype_used_in_interaction_focal_species;
		int phenotype_used_in_interaction_alternative_species;

		float mean_number_of_others_sampled;
		float deme_specific_interaction_effects_on_fecundity;
		float deme_specific_interaction_effects_on_survivorship;
		
		int FECUNDITY_PHENOTYPE_INDEX;
		int MORTALITY_PHENOTYPE_INDEX;

		void prepare_interactions(int alt_species_index);
		void calculate_cumulative_deme_sizes();
		void interact(int target_species);
		void survive();
	};

struct survivorship_functor
{
	/* 
		Elements in the tuple.

		----------------------
		0: cumulative effects of interactions
		1: phenotype governing survivorship
	
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float fitness_effect = 1.0 + thrust::get<0>(t);
		thrust::get<1>(t) = fitness_effect*thrust::get<1>(t);
		}
};
#endif
