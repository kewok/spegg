#include <util/amplify.h>
#include <util/which_function.h>


/* If you have a vector of values

A = [A0, A1, A2, A3, ..., An]

and each element gets repeated


B = [B0, B1, B2, B3, ..., Bn] 

times, produce a third vector C that repeats A0 B0 times, A1 B1 times, etc...

e.g., 

A = [0, 1, 2, 3, 4]

B = [5, 2, 5, 0, 1]


then 

C = [0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 4]


*/


/********************************

This function amplify takes two vectors, amplify_counts and values_to_amplify, and stores the results in a vector amplified_values where values_to_amplify[i] is repeated amplify_counts[i] times.

E.g., 

thrust::device_vector<int> amplify_counts(5);
thrust::device_vector<int> values_to_amplify(5);
thrust::device_vector<int> amplified_values(1);

for (int i=0; i < 5; i++)
	{
	values_to_amplify[i] = i*2;
	amplify_counts[i] = 4 + i;
	}

amplify(values_to_amplify, amplify_counts, amplified_values);

# amplified_values = [0, 0, 0, 0, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
	
********************************/

void amplify(thrust::device_vector<int> &values_to_amplify,
	     thrust::device_vector<int> &amplify_counts,
	     thrust::device_vector<int> &amplified_values)
	{
	thrust::device_vector<int> values_to_amplify_copy(values_to_amplify.size());
	thrust::device_vector<int> amplify_counts_copy(amplify_counts.size());

	/* these next 4 steps are needed to make sure that values with 0 amplifications are not in the final version */
	thrust::device_vector<int> viable_values;
	which_greater_than(amplify_counts, viable_values, 0);

	thrust::gather(viable_values.begin(), viable_values.end(), values_to_amplify.begin(), values_to_amplify_copy.begin());

	thrust::gather(viable_values.begin(), viable_values.end(), amplify_counts.begin(), amplify_counts_copy.begin());

	int number_of_values_to_amplify = values_to_amplify_copy.size();
	int total_number_of_amplified_values = thrust::reduce(amplify_counts_copy.begin(), amplify_counts_copy.end());

	amplified_values.resize(total_number_of_amplified_values);

	thrust::device_vector<int> amplify_counts_offsets( number_of_values_to_amplify );
	thrust::inclusive_scan(amplify_counts_copy.begin(), amplify_counts_copy.end(), amplify_counts_offsets.begin()); 

	thrust::device_vector<int> temporary_sequence ( total_number_of_amplified_values );
	thrust::sequence( temporary_sequence.begin(), temporary_sequence.end());

	thrust::upper_bound(amplify_counts_offsets.begin(), amplify_counts_offsets.end(), temporary_sequence.begin(), temporary_sequence.end(), amplified_values.begin());

	thrust::gather(amplified_values.begin(), amplified_values.end(), values_to_amplify_copy.begin() , amplified_values.begin());
	}

// Note there should actually be a step that checks to remove values where B=0 before proceeding.
void amplify_float(thrust::device_vector<float> &values_to_amplify,
		   thrust::device_vector<int> &amplify_counts,
	           thrust::device_vector<float> &amplified_values)
	{
	int number_of_values_to_amplify = values_to_amplify.size();
	int total_number_of_amplified_values = thrust::reduce(amplify_counts.begin(), amplify_counts.end());

	amplified_values.resize(total_number_of_amplified_values);

	thrust::device_vector<int> amplify_counts_offsets( number_of_values_to_amplify );
	thrust::inclusive_scan(amplify_counts.begin(), amplify_counts.end(), amplify_counts_offsets.begin()); 

	thrust::device_vector<float> temporary_sequence ( total_number_of_amplified_values );
	thrust::sequence( temporary_sequence.begin(), temporary_sequence.end());

	thrust::upper_bound(amplify_counts_offsets.begin(), amplify_counts_offsets.end(), temporary_sequence.begin(), temporary_sequence.end(), amplified_values.begin());

	thrust::gather(amplified_values.begin(), amplified_values.end(), values_to_amplify.begin() , amplified_values.begin());
	}

void amplify_sequence(thrust::device_vector<int> &amplify_counts,
		      int number_of_elements_in_sequence,
		      thrust::device_vector<int> &amplified_values)
	{
	thrust::device_vector<int> values_to_amplify(number_of_elements_in_sequence);
	thrust::sequence(values_to_amplify.begin(), values_to_amplify.end());

	int number_of_values_to_amplify = values_to_amplify.size();
	int total_number_of_amplified_values = thrust::reduce(amplify_counts.begin(), amplify_counts.end());

	amplified_values.resize(total_number_of_amplified_values);

	thrust::device_vector<int> amplify_counts_offsets( number_of_values_to_amplify );
	thrust::inclusive_scan(amplify_counts.begin(), amplify_counts.end(), amplify_counts_offsets.begin()); 

	thrust::device_vector<int> temporary_sequence ( total_number_of_amplified_values );
	thrust::sequence( temporary_sequence.begin(), temporary_sequence.end());

	thrust::upper_bound(amplify_counts_offsets.begin(), amplify_counts_offsets.end(), temporary_sequence.begin(), temporary_sequence.end(), amplified_values.begin());

	thrust::gather(amplified_values.begin(), amplified_values.end(), values_to_amplify.begin() , amplified_values.begin());
	}
		
