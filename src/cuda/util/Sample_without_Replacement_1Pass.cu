#include <util/Sample_without_Replacement_1Pass.h>
#include <util/amplify.h>
#include <util/thrust_functors.h>
#include <util/Sample_With_Replacement.h>
#include <util/remove_duplicate_pairs.h>

void Sample_without_Replacement_1Pass::sample()
	{
	// Reset sampling_input to allow sampling with replacement for one round:
	sampling_input->sampling_scheme = "SAMPLE_WITH_REPLACEMENT";

	class SamplingEvent *sample_others;
	sample_others = sample_others->create_SamplingEvent(sampling_input, gen);
	sample_others->sample();

	sampling_individuals_indices.resize(sample_others->sampling_individuals_indices.size());
	sampled_individuals_indices.resize(sample_others->sampled_individuals_indices.size());

	thrust::copy(sample_others->sampling_individuals_indices.begin(), sample_others->sampling_individuals_indices.end(), sampling_individuals_indices.begin());

	thrust::copy(sample_others->sampled_individuals_indices.begin(), sample_others->sampled_individuals_indices.end(), sampled_individuals_indices.begin());

	// Remove duplicates
	remove_duplicate_pairs(sampling_individuals_indices, sampled_individuals_indices);
	}

