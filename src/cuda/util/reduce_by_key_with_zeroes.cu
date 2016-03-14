#include <util/reduce_by_key_with_zeroes.h>

/* This function extends the output from thrust::reduce_by_key().


thrust::reduce_by_key() will not return a "zero" for values for which a key does not exist. For example,

	max_possible_keys = 4

	thrust::device_vector<int> key(5);
	thrust::device_vector<int> val(5);
	
	key[0] = 1;
	key[1] = 1;
	key[2] = 1;
	key[3] = 2;
	key[4] = 2;

	for (int i=0; i < 5; i++)
		val[i] = i;

	thrust::device_vector<int> key_out(max_possible_keys);
	thrust::device_vector<int> output(max_possible_keys);

	thrust::reduce_by_key(key.begin(), key.end(), val.begin(), key_out.begin(), output.begin());	

returns, in the vector output = [3, 7] and key_out=[1,2].

However, it does not give us output_complete = [0, 3, 7, 0] for keys =[0,3]. This function returns output_complete.

*/

void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<int> &values, thrust::device_vector<int> &values_output, int number_of_values, int number_of_possible_keys)
	{	
	// Figure out the number m of unique elements and store the values of each unique element in key_copy
	// e.g., if keys are [1,1,2,4,4,4] then key_copy should be [1,2,4]
	thrust::device_vector<int> key_copy(number_of_possible_keys);
	key_copy.erase(thrust::unique_copy(key_vector.begin(), key_vector.begin() + number_of_values, key_copy.begin()), key_copy.end());

	int number_of_keys_available = thrust::distance(key_copy.begin(), key_copy.end());

	thrust::device_vector<int> key_outputs(number_of_keys_available);
	thrust::device_vector<int> reduction_results(number_of_keys_available);

	thrust::reduce_by_key(key_vector.begin(), key_vector.begin() + number_of_values, values.begin(), key_outputs.begin(), reduction_results.begin()); 

	values_output.resize(number_of_possible_keys);

	thrust::fill(values_output.begin(), values_output.end(), 0);

	thrust::scatter(reduction_results.begin(), reduction_results.end(), key_outputs.begin(), values_output.begin());
	}


void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<float> &values, thrust::device_vector<float> &values_output, int number_of_values, int number_of_possible_keys)
	{	
	// Figure out the number m of unique elements and store the values of each unique element in key_copy
	// e.g., if keys are [1,1,2,4,4,4] then key_copy should be [1,2,4]

	thrust::device_vector<int> key_copy(number_of_possible_keys);
	key_copy.erase(thrust::unique_copy(key_vector.begin(), key_vector.begin() + number_of_values, key_copy.begin()), key_copy.end());
	int number_of_keys_available = thrust::distance(key_copy.begin(), key_copy.end());

	thrust::device_vector<int> key_outputs(number_of_keys_available);
	thrust::device_vector<float> reduction_results(number_of_keys_available);

	thrust::reduce_by_key(key_vector.begin(), key_vector.begin() + number_of_values, values.begin(), key_outputs.begin(), reduction_results.begin()); 

	values_output.resize(number_of_possible_keys);

	thrust::fill(values_output.begin(), values_output.end(), 0);

	thrust::scatter(reduction_results.begin(), reduction_results.end(), key_outputs.begin(), values_output.begin());
	}

