#include <util/Sample_without_Replacement_Test.h>
#include <util/thrust_functors.h>
#include <util/amplify.h>

#include <thrust/unique.h>

// This class draws a single sample, sampled without replacement, for each individual potentially conducting the sampling.

void Sample_without_Replacement_Test::sample()
	{
	 sampled_individuals_indices.resize(number_of_individuals_subject_to_sampling);
	thrust::copy(sampling_input->list_of_individuals_potentially_subject_to_sampling.begin(), sampling_input->list_of_individuals_potentially_subject_to_sampling.end(),  sampled_individuals_indices.begin());

	shuffle_sampled_individuals_at_random();

	amplify(sampling_input->list_of_individuals_potentially_conducting_sampling, sampling_input->number_of_other_individuals_sampled, sampling_individuals_indices);

	sampled_individuals_indices.erase(sampled_individuals_indices.begin() + sampling_individuals_indices.size(), sampled_individuals_indices.end());
	}

void Sample_without_Replacement_Test::shuffle_sampled_individuals_at_random()
	{
	unique_uniform_rvs.resize(number_of_individuals_subject_to_sampling);

	int max_number_of_duplicates = 1;

	while(max_number_of_duplicates > 0)
		{
		// This algorithm can be interminable if number_of_individuals_subject_to_sampling is huge (not clear what magic number is; testing seems to suggest it still handles about 70 million random variates fine, which is at the memory limit of the GPU)
		double *rand_ptr = raw_pointer_cast(&unique_uniform_rvs[0]);
		curandGenerateUniformDouble(gen, rand_ptr, number_of_individuals_subject_to_sampling); 

		thrust::device_vector<double> sorted_rand(number_of_individuals_subject_to_sampling);
		thrust::copy(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), sorted_rand.begin());
		thrust::sort(sorted_rand.begin(), sorted_rand.end());
		// If there are duplicates, this line of code should remove them from the sorted_rand vector.
		sorted_rand.erase(thrust::unique(sorted_rand.begin(), sorted_rand.end()), sorted_rand.end());
		// Consequently, the sorted_rand vector will be shorter than the original uniform RV vector, necessitating a new set of uniform RVs be drawn:
		max_number_of_duplicates = unique_uniform_rvs.size() - sorted_rand.size();
		
		if (max_number_of_duplicates == 0)
			{
			thrust::stable_sort_by_key(unique_uniform_rvs.begin(), unique_uniform_rvs.end(),  sampled_individuals_indices.begin());
			}
		}
	}
