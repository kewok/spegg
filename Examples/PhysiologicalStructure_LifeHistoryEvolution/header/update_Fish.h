#ifndef UPDATE_FISH_H
#define UPDATE_FISH_H

#include <species/update/updatebehavior.h>
#include "Fish_Habitat.h"

#include <thrust/sequence.h>
#include <thrust/transform.h>

class update_Fish : public UpdateBehavior
	{
	inds_stochastic *species;
	Fish_Habitat *habitat;

	// Constructor
	public:
	update_Fish(inds_stochastic *species, environment *habitat) 
	 	{
		this->species = species;

		this->MORTALITY_PHENOTYPE = species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE"];
		this->FECUNDITY_PHENOTYPE = species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE"];
		this->IRREVERSIBLE_MASS_PHENOTYPE = species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];
		this->REVERSIBLE_MASS_PHENOTYPE = species->demeParameters->species_specific_values["REVERSIBLE_MASS_PHENOTYPE"];
		this->RESOURCE_LIMITATION_PHENOTYPE = species->demeParameters->species_specific_values["RESOURCE_LIMITATION_PHENOTYPE"];
		this->EGGSIZE_PHENOTYPE =  species->demeParameters->species_specific_values["EGGSIZE_PHENOTYPE"];

		// Copy the constants 
		this->size = species->size;
		this->Number_of_Demes = species->Num_Demes;

		// Nifty hack: Downcast your generic habitat object to a pointer to fish class to access its members
		this->habitat = static_cast<Fish_Habitat*> (habitat);
		this->intra_annual_time_steps = this->habitat->intra_annual_time_steps;

		// Prepare pointers to thrust vectors:
		prepare_growth_constants_pointers();
		prepare_survivorship_constants_pointers();
		}

		void update();

	protected:
		// variable names
		int MORTALITY_PHENOTYPE;
		int FECUNDITY_PHENOTYPE;
		int IRREVERSIBLE_MASS_PHENOTYPE;
		int REVERSIBLE_MASS_PHENOTYPE;
		int RESOURCE_LIMITATION_PHENOTYPE;
		int EGGSIZE_PHENOTYPE;

		int size;
		int Number_of_Demes;
		int intra_annual_time_steps;
		void consume();
		void survive();
		void update_vital_rates();
		void update_fecundity();

		void prepare_growth_constants_pointers();
		void prepare_survivorship_constants_pointers();

		// Pointers to thrust vectors used in survivorship calculation:
		float *effect_of_starvation_ptr;
		float *size_dependent_mortality_constant_ptr;
		float *size_dependent_mortality_coefficient_ptr;

		// Pointers to thrust vectors used in somatic growth calculation:
		float *consumption_allometric_scalar_ptr;
		float *gamma_ptr;
		float *ontogenetic_niche_shift_constant_ptr;
		float *ontogenetic_niche_shift_coefficient_ptr;
		float *length_weight_conversion_coefficient_ptr;
		float *length_weight_conversion_exponent_ptr;

		float *resource_1_maximum_ptr;
		float *resource_2_maximum_ptr;
		float *handling_time_resource1_ptr;
		float *handling_time_resource2_ptr;

		float *functional_response_scalar_resource1_ptr;
		float *functional_response_scalar_resource2_ptr;
		float *alpha_g_ptr;

		float *mature_maximum_condition_ptr;
		float *juvenile_maximum_condition_ptr;

		float *M_sizes_at_maturity_ptr;
		float *F_sizes_at_maturity_ptr;
	};
#endif
