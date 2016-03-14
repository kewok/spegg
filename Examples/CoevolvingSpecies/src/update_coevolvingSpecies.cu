#include "update_coevolvingSpecies.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void update_coevolvingSpecies::update()
	{
	for (int i=0; i < species[FOCAL_SPECIES_INDEX]->demeParameters->interacting_species.size(); i++)
		{
		prepare_interactions(i);
		interact(i);
		}
	}


// This nomenclature is confusing and should be fixed; target_species refers to the index AMONG interacting species, whilst ALTERNATIVE_SPECIES_INDEX refers to the index of target_species among all species.
void update_coevolvingSpecies::prepare_interactions(int target_species)
	{
	ALTERNATIVE_SPECIES_INDEX = species[FOCAL_SPECIES_INDEX]->demeParameters->interacting_species[target_species];

	calculate_cumulative_deme_sizes();

	phenotype_used_in_interaction_focal_species = species[FOCAL_SPECIES_INDEX]->demeParameters->interaction_phenotype_indices[ALTERNATIVE_SPECIES_INDEX];

	phenotype_used_in_interaction_alternative_species = species[ALTERNATIVE_SPECIES_INDEX]->demeParameters->interaction_phenotype_indices[FOCAL_SPECIES_INDEX];

	mean_number_of_others_sampled = *(species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("mean_number_of_others_sampled"));
	deme_specific_interaction_effects_on_fecundity = species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_fecundity[target_species][0];
	deme_specific_interaction_effects_on_survivorship = species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_survivorship[target_species][0];
	}

void update_coevolvingSpecies::calculate_cumulative_deme_sizes()
	{
	cumulative_deme_sizes.resize(species[FOCAL_SPECIES_INDEX]->Num_Demes);
	thrust::inclusive_scan(species[ALTERNATIVE_SPECIES_INDEX]->deme_sizes.begin(), species[ALTERNATIVE_SPECIES_INDEX]->deme_sizes.begin() + species[ALTERNATIVE_SPECIES_INDEX]->Num_Demes, cumulative_deme_sizes.begin());
	}

void update_coevolvingSpecies::interact(int target_species)
	{
	int deme_offset = 0;
	
	float fecundity_trait_trade_off = *(species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("effect_of_interaction_phenotype_on_fecundity"));
	float effect_of_interaction_on_fecundity = 1;
	float effect_of_interaction_on_survivorship = 1;

	int min_ind_index = 0;
	int max_ind_index = cumulative_deme_sizes[0];

	// Miscellaneous helper variables:
	int number_of_others_sampled, alternative_individual;
	float use_rv, weighted_effect_of_alternative_individual, old_fecundity, old_survivorship;

	for (int i=0; i < species[FOCAL_SPECIES_INDEX]->size; i++)
		{
		if (species[FOCAL_SPECIES_INDEX]->status[i] > 0)// If the individual is alive; this should be universally true; so no indentation for clarity.
		{
		// identify the deme_specific values; optimized to minimize memory access:
		if (species[FOCAL_SPECIES_INDEX]->deme[i] != deme_offset)
			{
			deme_offset = species[FOCAL_SPECIES_INDEX]->deme[i];
			mean_number_of_others_sampled = *(species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("mean_number_of_others_sampled") + deme_offset);

			deme_specific_interaction_effects_on_fecundity = species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_fecundity[target_species][deme_offset];
			deme_specific_interaction_effects_on_survivorship = species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_survivorship[target_species][deme_offset];

			fecundity_trait_trade_off = *(species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("effect_of_interaction_phenotype_on_fecundity") + deme_offset);
	
			min_ind_index = cumulative_deme_sizes[deme_offset - 1];
			max_ind_index = cumulative_deme_sizes[deme_offset];
			}

		float focal_individuals_phenotype = species[FOCAL_SPECIES_INDEX]->phenotype[phenotype_used_in_interaction_focal_species][i];

		if (min_ind_index != max_ind_index) // provided there is at least one other individual that is alive:
			{
			effect_of_interaction_on_fecundity = 1;
			effect_of_interaction_on_survivorship = 1;

			if (mean_number_of_others_sampled > MAX_SAMPLED)
				{
				// Use the discrete normal approximation to the poisson
				number_of_others_sampled = (int) round((double) mean_number_of_others_sampled + gsl_ran_gaussian(species[FOCAL_SPECIES_INDEX]->gen, sqrt(mean_number_of_others_sampled)));
				if (number_of_others_sampled < 0)
					number_of_others_sampled = 0;
				}
			else
				{
				number_of_others_sampled = gsl_ran_poisson(species[FOCAL_SPECIES_INDEX]->gen, mean_number_of_others_sampled);
				}

			for (int k=0; k < number_of_others_sampled; k++)
				{
				use_rv = gsl_rng_uniform(species[FOCAL_SPECIES_INDEX]->gen);
				alternative_individual = draw_ind(min_ind_index, max_ind_index, use_rv);
				weighted_effect_of_alternative_individual = 1;

				// Discount effect if the individual from the alternative species is dead, or if the individual sampled themselves. Note that alternative suitability criteria can be used here; e.g., if the alternative individual is a male, if the alternative individual is yourself, etc...
				if ( (species[ALTERNATIVE_SPECIES_INDEX]->status[alternative_individual] == 0) || ((FOCAL_SPECIES_INDEX == ALTERNATIVE_SPECIES_INDEX) && (alternative_individual == i) ) )
					{
					weighted_effect_of_alternative_individual = 0;	
					}

// For heterospecifics, the interaction strengths are mediated by trait values focal_phen1 and alt_phen1
				if ((FOCAL_SPECIES_INDEX != ALTERNATIVE_SPECIES_INDEX))
					{
					effect_of_interaction_on_fecundity = calculate_effect_of_interaction(effect_of_interaction_on_fecundity, weighted_effect_of_alternative_individual, deme_specific_interaction_effects_on_fecundity, focal_individuals_phenotype, species[ALTERNATIVE_SPECIES_INDEX]->phenotype[phenotype_used_in_interaction_alternative_species][alternative_individual]);
					effect_of_interaction_on_survivorship = calculate_effect_of_interaction(effect_of_interaction_on_survivorship, weighted_effect_of_alternative_individual, deme_specific_interaction_effects_on_survivorship, focal_individuals_phenotype, species[ALTERNATIVE_SPECIES_INDEX]->phenotype[phenotype_used_in_interaction_alternative_species][alternative_individual]);
					}
				else
					{
					effect_of_interaction_on_fecundity = calculate_effect_of_interaction_noPhen(effect_of_interaction_on_fecundity, weighted_effect_of_alternative_individual, deme_specific_interaction_effects_on_fecundity);
					effect_of_interaction_on_survivorship = calculate_effect_of_interaction_noPhen(effect_of_interaction_on_survivorship, weighted_effect_of_alternative_individual, deme_specific_interaction_effects_on_survivorship);
					}
				}
			
			}
		old_fecundity = species[FOCAL_SPECIES_INDEX]->phenotype[FECUNDITY_PHENOTYPE_INDEX][i];
		species[FOCAL_SPECIES_INDEX]->phenotype[FECUNDITY_PHENOTYPE_INDEX][i] = update_fecundity_phenotype(effect_of_interaction_on_fecundity, old_fecundity, focal_individuals_phenotype, fecundity_trait_trade_off);

		old_survivorship = species[FOCAL_SPECIES_INDEX]->phenotype[MORTALITY_PHENOTYPE_INDEX][i];
		species[FOCAL_SPECIES_INDEX]->phenotype[MORTALITY_PHENOTYPE_INDEX][i] = update_survivorship_phenotype(effect_of_interaction_on_survivorship, old_survivorship);
		}
		}
	}
