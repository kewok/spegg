#include <math/histogram.h>
#include <util/Sample_without_Replacement_byDeme.h>

void Sample_without_Replacement_byDeme::sample()
	{
	int max_number_of_duplicates = 1;

	/* For each of the individuals conducting the sampling, identify the index of the shuffled individual to which they would be assigned */
	int *cumulative_demewise_sum_species_conducting_sampling = raw_pointer_cast(&cumulative_sampling_individuals_by_deme[0]);
	int *cumulative_demewise_sum_species_subject_to_sampling = raw_pointer_cast(&cumulative_sampleable_individuals_by_deme[0]);

	thrust::device_vector<int> sampling_individuals(sampling_individuals_indices.size());
	thrust::sequence(sampling_individuals.begin(), sampling_individuals.end());
	// Identify the index of the individual to be sampled, stratified by deme
	specify_index_to_sample idx_to_sample(cumulative_demewise_sum_species_conducting_sampling, cumulative_demewise_sum_species_subject_to_sampling);

	thrust::for_each(thrust::make_zip_iterator(
			 	thrust::make_tuple(sampling_individuals.begin(), sampling_input->deme_affiliation_of_sampling_individuals.begin(), index_to_sample.begin())), 
			 thrust::make_zip_iterator(
			 	thrust::make_tuple(sampling_individuals.end(), sampling_input->deme_affiliation_of_sampling_individuals.end(), index_to_sample.end())), idx_to_sample);

	cudaThreadSynchronize();

	while(max_number_of_duplicates > 0)
		{
		// Note this routine suffers from the same interminable problem outlined in Sample_without_Replacement_Test if there are identical uniform rvs. Test simulations suggest the scenarios that will cause you troubles is unlikely to arise frequently in practice with 10^7 double variates
		double *rand_ptr = raw_pointer_cast(&unique_uniform_rvs[0]);
		curandGenerateUniformDouble(gen, rand_ptr, sampling_input->list_of_individuals_potentially_subject_to_sampling.size());

		thrust::device_vector<double> temp_demes(demes_of_individuals_subject_to_sampling.size());

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(demes_of_individuals_subject_to_sampling.begin(), temp_demes.begin())),
				 thrust::make_zip_iterator(thrust::make_tuple(demes_of_individuals_subject_to_sampling.end(), temp_demes.end())), 
				 int_to_double());

		thrust::transform(temp_demes.begin(), temp_demes.end(), unique_uniform_rvs.begin(), unique_uniform_rvs.begin(), thrust::plus<double>());
		cudaThreadSynchronize();
		/* Test for any duplicates */
		thrust::device_vector<double> sorted_rand(sampling_input->list_of_individuals_potentially_subject_to_sampling.size());

		thrust::copy(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), sorted_rand.begin());
		thrust::sort(sorted_rand.begin(), sorted_rand.end());
		// If there are duplicates, this line of code should remove them from the sorted_rand vector.
		sorted_rand.erase(thrust::unique(sorted_rand.begin(), sorted_rand.end()), sorted_rand.end());
		// Consequently, the sorted_rand vector will be shorter than the original uniform RV vector, necessitating a new set of uniform RVs be drawn:
		max_number_of_duplicates = unique_uniform_rvs.size() - sorted_rand.size();	
		// Provided there were no duplicates
		if (max_number_of_duplicates == 0)
			{
			thrust::device_vector<int> temp_sampled_individuals_indices(sampling_input->list_of_individuals_potentially_subject_to_sampling.size());
			thrust::copy(sampling_input-> list_of_individuals_potentially_subject_to_sampling.begin(), sampling_input-> list_of_individuals_potentially_subject_to_sampling.end(), temp_sampled_individuals_indices.begin());
			thrust::stable_sort_by_key(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), temp_sampled_individuals_indices.begin());

			thrust::gather(index_to_sample.begin(), index_to_sample.end(), temp_sampled_individuals_indices.begin(),  sampled_individuals_indices.begin());
			sampled_individuals_indices.erase(sampled_individuals_indices.begin() + sampling_individuals_indices.size(), sampled_individuals_indices.end());
			}
		}
	}

void Sample_without_Replacement_byDeme::setup_demes(thrust::device_vector<int> &demes_of_individuals_sampled)
	{
	thrust::copy(demes_of_individuals_sampled.begin(), demes_of_individuals_sampled.end(), demes_of_individuals_subject_to_sampling.begin());
	calculate_histogram(demes_of_individuals_sampled, sampling_input->sampleable_individuals_per_deme, sampling_input->Num_Demes);

	thrust::exclusive_scan(sampling_input->sampleable_individuals_per_deme.begin(), sampling_input->sampleable_individuals_per_deme.end(), cumulative_sampleable_individuals_by_deme.begin());
	
	calculate_histogram(sampling_input->deme_affiliation_of_sampling_individuals, number_of_sampling_individuals_by_deme, sampling_input->Num_Demes);
	thrust::exclusive_scan(number_of_sampling_individuals_by_deme.begin(),  number_of_sampling_individuals_by_deme.end(), cumulative_sampling_individuals_by_deme.begin());
	}
