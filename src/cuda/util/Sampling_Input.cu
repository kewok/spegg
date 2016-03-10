#include "sampling_input.h"
#include "inds.h"
//#include "species_specific_sampling_input.h"
#include "histogram.h"

// The operations in prepare_sampling_input() are also shared across species and contexts
#if 0
void SamplingInput::prepare_sampling_input()
	{
	specify_individuals_conducting_sampling();
	specify_individuals_being_sampled();
	set_parameters();
	determine_number_of_individuals_sampled();
	}


// Methods to determine how many individuals each individual samples
void SamplingInput::determine_number_of_individuals_to_be_sampled_poisson() 
	{

	int number_of_individuals_doing_the_sampling = list_of_individuals_potentially_conducting_sampling.size();

	// For each individual, based on their deme, identify the expected number of other individuals they are likely to sample
	thrust::device_vector<float> mean_numbers_of_others_sampled( number_of_individuals_doing_the_sampling );
	
	thrust::gather(population_affiliation_of_sampling_individuals.begin(), 
		       population_affiliation_of_sampling_individuals.end(), 
		       species[focal_species_ID]->subpopParameters->get_vector_ptr("mean_number_of_others_sampled"),
		       mean_numbers_of_others_sampled.begin());

	// draw a poisson, where mean_numbers_of_others_sampled is the mean.
	draw_poisson_different_parameters(number_of_individuals_doing_the_sampling, 
					  mean_numbers_of_others_sampled, 
					  number_of_other_individuals_sampled,
					  species[focal_species_ID]->gen);

	// Make sure you sample no more than there are individuals in the other species deme by:
		// 1. Figuring out the demes of individuals subject to being sampled
	thrust::device_vector<int> maximum_sampleable_individuals( number_of_individuals_doing_the_sampling );

	thrust::device_vector<int> target_species_deme_values( list_of_individuals_potentially_subject_to_sampling.size() );

	thrust::gather(list_of_individuals_potentially_subject_to_sampling.begin(),
		       list_of_individuals_potentially_subject_to_sampling.end(),
		       species[target_species_ID]->pop.begin(),
		       target_species_deme_values.begin());
	
		// 2. Count the number of times each deme is represented
	thrust::device_vector<int> possible_deme_values( species[target_species_ID]->Num_Subpopulations );
	thrust::sequence(possible_deme_values.begin(), possible_deme_values.end());

	calculate_histogram(target_species_deme_values, sampleable_individuals_per_deme, species[target_species_ID]->Num_Subpopulations);

	thrust::gather(population_affiliation_of_sampling_individuals.begin(),
		       population_affiliation_of_sampling_individuals.end(),
		       sampleable_individuals_per_deme.begin(),
		       maximum_sampleable_individuals.begin());

		// 3. Apply correction for small deme size: cannot sample more individuals than there are individuals in the deme.
	thrust::transform( number_of_other_individuals_sampled.begin(),  number_of_other_individuals_sampled.end(), maximum_sampleable_individuals.begin(),  number_of_other_individuals_sampled.begin(), thrust::minimum<int>());
	}


// Method to set the same number of others sampled for all individuals in a deme
void SamplingInput::determine_number_of_individuals_to_be_sampled_fixed() 
	{
	// Assume all individuals sample the same number of individuals in the other species.
	thrust::device_vector<int> maximum_sampleable_individuals( number_of_individuals_in_sampling_species );
	thrust::device_vector<int> target_species_deme_values( list_of_individuals_potentially_subject_to_sampling.size() );

	thrust::gather(list_of_individuals_potentially_subject_to_sampling.begin(),
		       list_of_individuals_potentially_subject_to_sampling.end(),
		       species[target_species_ID]->pop.begin(),
		       target_species_deme_values.begin());
	
	thrust::device_vector<int> possible_deme_values( species[target_species_ID]->Num_Subpopulations );
	thrust::sequence(possible_deme_values.begin(), possible_deme_values.end());

	calculate_histogram(target_species_deme_values, sampleable_individuals_per_deme, species[target_species_ID]->Num_Subpopulations);

	thrust::gather(population_affiliation_of_sampling_individuals.begin(),
		       population_affiliation_of_sampling_individuals.end(),
		       sampleable_individuals_per_deme.begin(),
		       maximum_sampleable_individuals.begin());

	// specify the number of other individuals sampled, which we assume is set deterministically
	thrust::gather(population_affiliation_of_sampling_individuals.begin(), 
		       population_affiliation_of_sampling_individuals.end(), 
		       species[focal_species_ID]->subpopParameters->get_vector_ptr("mean_number_of_others_sampled"),
		       number_of_other_individuals_sampled.begin());
	

	// Apply correction for small deme size: cannot sample more individuals than there are individuals in the deme.
	thrust::transform( number_of_other_individuals_sampled.begin(),  number_of_other_individuals_sampled.end(), maximum_sampleable_individuals.begin(),  number_of_other_individuals_sampled.begin(), thrust::minimum<int>());
	}
#endif
