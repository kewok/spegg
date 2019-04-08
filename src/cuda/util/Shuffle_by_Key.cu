#include <util/shuffle_by_key.h>

void shuffle_by_key(thrust::device_vector<int> &keys, thrust::device_vector<int> &values, curandGenerator_t gen, bool shuffle_keys)
	{
	int max_number_of_duplicates = 1;
	thrust::device_vector<double> unique_uniform_rvs(keys.size());

	while(max_number_of_duplicates > 0)
		{
		// Note this routine suffers from the same interminable problem outlined in Sample_without_Replacement_Test if there are identical uniform rvs. Test simulations suggest the scenarios that will cause you troubles is unlikely to arise frequently in practice with 10^7 double variates
		double *rand_ptr = raw_pointer_cast(&unique_uniform_rvs[0]);
		curandGenerateUniformDouble(gen, rand_ptr, keys.size());

		thrust::device_vector<double> temp_keys(keys.size());

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(keys.begin(), temp_keys.begin())),
				 thrust::make_zip_iterator(thrust::make_tuple(keys.end(), temp_keys.end())), 
				 int_to_double());

		thrust::transform(temp_keys.begin(), temp_keys.end(), unique_uniform_rvs.begin(), unique_uniform_rvs.begin(), thrust::plus<double>());
		cudaThreadSynchronize();
		/* Test for any duplicates */
		thrust::device_vector<double> sorted_rand(keys.size());

		thrust::copy(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), sorted_rand.begin());
		thrust::sort(sorted_rand.begin(), sorted_rand.end());
		// If there are duplicates, this line of code should remove them from the sorted_rand vector.
		sorted_rand.erase(thrust::unique(sorted_rand.begin(), sorted_rand.end()), sorted_rand.end());
		// Consequently, the sorted_rand vector will be shorter than the original uniform RV vector, necessitating a new set of uniform RVs be drawn:
		max_number_of_duplicates = unique_uniform_rvs.size() - sorted_rand.size();	
		// Provided there were no duplicates
		if (max_number_of_duplicates == 0)
			{
			if (shuffle_keys)
				{
				thrust::device_vector<double> unique_uniform_copy(keys.size());
				thrust::copy(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), unique_uniform_copy.begin());
				thrust::stable_sort_by_key(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), values.begin());
				thrust::stable_sort_by_key(unique_uniform_copy.begin(), unique_uniform_copy.end(), keys.begin());
				}
			else
				{
				thrust::stable_sort_by_key(unique_uniform_rvs.begin(), unique_uniform_rvs.end(), values.begin());
				}
			}
		}
	}
