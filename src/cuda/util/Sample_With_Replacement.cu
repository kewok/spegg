#include <util/Sample_With_Replacement.h>
#include <util/amplify.h>
#include <util/thrust_functors.h>

void Sample_With_Replacement::sample()
	{
	/*************
	//
	// 1. For each sampling event, generate a vector "individuals_sampling" that repeats each index i in (0,1,2,... , |x|) n_i times.
	//
	**************/

	amplify(sampling_input->list_of_individuals_potentially_conducting_sampling, sampling_input->number_of_other_individuals_sampled, sampling_individuals_indices);

	/*************
	//
	// 2. Perform the stratified sampling by deme
	//
	**************/

		/***********
		//
		// 2.a Generate a vector u, where the jth element of u is the cumulative number of individuals subject to sampling in deme j
		//
		***********/

	// The individuals subject to sampling have to fall within the ranges of the total number of individuals in the sampler's deme. Also make sure Num_Demes is initialized; if is set to zero you will get a memory error.


	thrust::device_vector<int> cumulative_sampleable_individuals_by_deme( sampling_input->Num_Demes );

	thrust::inclusive_scan( sampling_input->sampleable_individuals_per_deme.begin(), sampling_input-> sampleable_individuals_per_deme.end() , cumulative_sampleable_individuals_by_deme.begin() );	

	int *cumulative_others_by_deme_ptr = raw_pointer_cast(&cumulative_sampleable_individuals_by_deme[0]);

	thrust::device_vector<int> deme_of_samplers;

	amplify(sampling_input->deme_affiliation_of_sampling_individuals, sampling_input->number_of_other_individuals_sampled, deme_of_samplers);

		/***********
		//
		// 2.b. For each individual i that is doing the sampling, draw a random float z[i] between 0 and 1
		//
		***********/

	thrust::device_vector<float> sampled_floats( sampling_individuals_indices.size() ); 
	float *curand_samples_ptr = raw_pointer_cast( &sampled_floats[0] );
	curandGenerateUniform(gen, curand_samples_ptr, sampling_individuals_indices.size());

		/***********
		//
		// 2.c. Let that random number fall between u[j-1], u[j] where j is the deme of the ith individual sampler.
		//
		***********/

	// sample total_others_sampled other individuals at random
	thrust::device_vector<int> individuals_sampled( sampling_individuals_indices.size() ); 

	// Use the functor reassign to change the random integer to fall within the range of individuals in the deme
	reassign_functor reassign(cumulative_others_by_deme_ptr);

	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(
					sampled_floats.begin(),
					deme_of_samplers.begin(),
					individuals_sampled.begin()
					)),
			  thrust::make_zip_iterator(thrust::make_tuple(
					sampled_floats.end(),
					deme_of_samplers.end(),
					individuals_sampled.end())),
				reassign
				);

	/*************
	//
	// 3. Map Z[i] to the vector w to produce a vector W, the list of individuals that get sampled for each sampling event  (the gather operation)
	//
	**************/

	// Specify the index of the sampled individual whose index in individuals_sampled corresponds to:	
	// Make sure you use the index from the overall deme to which the sampled individuals belong. 

	sampled_individuals_indices.resize( sampling_individuals_indices.size() );
	thrust::gather( individuals_sampled.begin(), individuals_sampled.end(), sampling_input-> list_of_individuals_potentially_subject_to_sampling.begin(),  sampled_individuals_indices.begin());
	}

