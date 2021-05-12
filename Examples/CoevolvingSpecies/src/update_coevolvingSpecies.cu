#include "update_coevolvingSpecies.h"

void update_coevolvingSpecies::update()
	{
	for (int i=0; i < species[FOCAL_SPECIES_INDEX]->demeParameters->interacting_species.size(); i++)
		{
		prepare_interactions(i);
		interact();
		}
	}

// This nomenclature is confusing and should be fixed; target_species refers to the index AMONG interacting species, whilst ALTERNATIVE_SPECIES_INDEX refers to the index of target_species among all species.


void update_coevolvingSpecies::prepare_interactions(int target_species)
	{	
	ALTERNATIVE_SPECIES_INDEX = species[FOCAL_SPECIES_INDEX]->demeParameters->interacting_species[target_species];

	calculate_cumulative_deme_sizes();

	phenotype_used_in_interaction_focal_species = species[FOCAL_SPECIES_INDEX]->demeParameters->interaction_phenotype_indices[ALTERNATIVE_SPECIES_INDEX];
	phenotype_used_in_interaction_alternative_species = species[ALTERNATIVE_SPECIES_INDEX]->demeParameters->interaction_phenotype_indices[FOCAL_SPECIES_INDEX];

	// cast the vector pointers
	focal_species_phenotype = raw_pointer_cast(&species[FOCAL_SPECIES_INDEX]->phenotype[phenotype_used_in_interaction_focal_species][0]);
	alt_species_phenotype = raw_pointer_cast(&species[ALTERNATIVE_SPECIES_INDEX]->phenotype[phenotype_used_in_interaction_alternative_species][0]);
	alt_species_status = raw_pointer_cast(&species[ALTERNATIVE_SPECIES_INDEX]->status[0]);
	cumulative_deme_sizes_ptr = raw_pointer_cast(&cumulative_deme_sizes[0]);
	
	indices.resize(size);
	thrust::sequence(indices.begin(), indices.end());
	
	mean_number_of_others_sampled.resize(size);
	interaction_effects_on_fecundity.resize(size);
	interaction_effects_on_survivorship.resize(size);
	
	// Fan out the mean number of individuals that each individual of FOCAL_SPECIES will encounter
	thrust::gather(species[FOCAL_SPECIES_INDEX]->deme.begin(), species[FOCAL_SPECIES_INDEX]->deme.begin() + size, species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("mean_number_of_others_sampled"), mean_number_of_others_sampled.begin());

	// Fan out the interaction effects on fecundity by deme
	thrust::gather(species[FOCAL_SPECIES_INDEX]->deme.begin(), species[FOCAL_SPECIES_INDEX]->deme.begin() + size, species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_fecundity[target_species].begin(), interaction_effects_on_fecundity.begin());

	// Fan out the interaction effects on survivorship by deme
	thrust::gather(species[FOCAL_SPECIES_INDEX]->deme.begin(), species[FOCAL_SPECIES_INDEX]->deme.begin() + size, species[FOCAL_SPECIES_INDEX]->demeParameters->deme_wise_interaction_effects_on_survivorship[target_species].begin(), interaction_effects_on_survivorship.begin());

	// Unit testing code:
/* 	for (int i=0; i < size; i++)
		std::cout << " check that things were correctly setup for interaction between " << FOCAL_SPECIES_INDEX << " " << ALTERNATIVE_SPECIES_INDEX << " " << species[FOCAL_SPECIES_INDEX]->deme[i] << " " << interaction_effects_on_survivorship[i] << " " << interaction_effects_on_fecundity[i] <<  std::endl;
*/
	}

void update_coevolvingSpecies::interact()
	{
	interaction_kernel upit(FOCAL_SPECIES_INDEX, ALTERNATIVE_SPECIES_INDEX, focal_species_phenotype, alt_species_phenotype, alt_species_status, cumulative_deme_sizes_ptr);

	/* randomize the seeds */
	thrust::device_vector<unsigned int> seed_vals(size);
	unsigned int *seed_ptr = raw_pointer_cast(&seed_vals[0]);
	curandGenerate(species[FOCAL_SPECIES_INDEX]->gen, seed_ptr, size);
	
	/* determine the cost of the phenotype - i.e., the trade-off strength */
	thrust::device_vector<float> cost_of_interaction_phenotype(size);
	thrust::gather(species[FOCAL_SPECIES_INDEX]->deme.begin(), species[FOCAL_SPECIES_INDEX]->deme.begin() + size, species[FOCAL_SPECIES_INDEX]->demeParameters->get_vector_ptr("effect_of_interaction_phenotype_on_fecundity"), cost_of_interaction_phenotype.begin());

	/* perform the simulation */
	thrust::for_each(thrust::make_zip_iterator(
					 thrust::make_tuple(indices.begin(), 								    species[FOCAL_SPECIES_INDEX]->deme.begin(), 
							    mean_number_of_others_sampled.begin(),
							    species[FOCAL_SPECIES_INDEX]->status.begin(), 
							    species[FOCAL_SPECIES_INDEX]->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),
							    species[FOCAL_SPECIES_INDEX]->phenotype[MORTALITY_PHENOTYPE_INDEX].begin(),
							    interaction_effects_on_fecundity.begin(), 
							    interaction_effects_on_survivorship.begin(),
							    cost_of_interaction_phenotype.begin(),
							    seed_vals.begin())),	
				 thrust::make_zip_iterator(
						thrust::make_tuple(indices.end(), 
							    species[FOCAL_SPECIES_INDEX]->deme.begin() + size,
							    mean_number_of_others_sampled.end(), 
							    species[FOCAL_SPECIES_INDEX]->status.begin() + size, 
							    species[FOCAL_SPECIES_INDEX]->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + size,
							    species[FOCAL_SPECIES_INDEX]->phenotype[MORTALITY_PHENOTYPE_INDEX].begin() + size,									   
							    interaction_effects_on_fecundity.end(),
							    interaction_effects_on_survivorship.end(),
							    cost_of_interaction_phenotype.end(),
							    seed_vals.end())),
				 upit);

	cudaDeviceSynchronize();
	}

void update_coevolvingSpecies::calculate_cumulative_deme_sizes()
	{
	cumulative_deme_sizes.resize(species[FOCAL_SPECIES_INDEX]->Num_Demes);
	thrust::inclusive_scan(species[ALTERNATIVE_SPECIES_INDEX]->deme_sizes.begin(), species[ALTERNATIVE_SPECIES_INDEX]->deme_sizes.begin() + species[FOCAL_SPECIES_INDEX]->Num_Demes, cumulative_deme_sizes.begin());
	}
