#include "mating_thrust_prob_table_demes.h"
#include "footimer2.h"
#include "thrust_functors.h"

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/gather.h>
#include <thrust/unique.h>

/*

Here is the idea with this implementation of thrust_prob_table

when we have keys [1,1,2,3,4,4,4] and probabilities [0.2,0.3,0.1,0.6,0.1,0.25,0.3], do an inclusive scan on the vector of probabilities, i.e., 
prob' = [0.20, 0.50, 0.60, 1.20, 1.30, 1.55, 1.85]
Do a histogram of [1, 1, 2, 3, 4, 4, 4] to get vector hist = [2, 1, 1, 3] (this is presumably done ahead of time, and the histogram is fed in as an argument to determine_key_offsets)
take an inclusive scan of hist, hist' = [2, 3, 4, 7]. These are your offsets.

Then for each uniform random variable u, if u is associated with key[i], u' = prob'[hist'[key[i]-1] ] + u*(prob'[hist'[key[i]] ]  - prob'[hist'[key[i]-1] ] )
Then draw based on u.
*/

// make sure the kid's uniform random number falls in the interval corresponding to their current population
struct adjust_randoms_functor
	{
	float *popRanges_in_cumulative;

	adjust_randoms_functor(float *pop_range) : popRanges_in_cumulative(pop_range)
	{};
	
	/*
		Elements in the tuple.
		----------------------
		0: the kid's random number
		1: the kid's population
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
			int this_kids_pop = thrust::get<1>(t);
			if (this_kids_pop != 0)
				{
				float interval = popRanges_in_cumulative[this_kids_pop] - popRanges_in_cumulative[this_kids_pop - 1];
				thrust::get<0>(t) = thrust::get<0>(t)*interval +  popRanges_in_cumulative[this_kids_pop-1];
				}

	 // For the case of the kids in population zero

			else
				{
				float interval = popRanges_in_cumulative[this_kids_pop];
				thrust::get<0>(t) = thrust::get<0>(t)*interval;
				}
		}
	};

void mating_ThrustProbTable_demes::determine_key_offsets(int number_of_key_types, thrust::device_vector<int> &key_histogram_vector )
	{
	thrust::device_vector<int> temp_offsets( number_of_key_types ); 
	thrust::inclusive_scan(key_histogram_vector.begin(), key_histogram_vector.end(), temp_offsets.begin());

	// Need to subtract one from the offsets (as offsets are based on counts, but we need to match index)
	key_offsets.resize( number_of_key_types );	
	thrust::transform(temp_offsets.begin(), temp_offsets.begin() + number_of_key_types, key_offsets.begin(), unary_minus<unsigned int>(1));
	}

void mating_ThrustProbTable_demes::adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end,
thrust::device_vector<int>::iterator inds_demes_begin, thrust::device_vector<int>::iterator inds_demes_end)
	{
	// Note the number of populations.

	int n = thrust::distance(key_offsets.begin(), key_offsets.end());

	// Figure out the bounds for each subpopulation in terms of their values in the cumulative probability table
	thrust::device_vector<float> bounds(n);

	thrust::gather(key_offsets.begin(), key_offsets.end(),cumulative_prob.begin(),bounds.begin());
	
	thrust::device_vector<int> temp(n);
	thrust::copy(key_offsets.begin(), key_offsets.end(), temp.begin());

	float *cumul_prob_bounds = raw_pointer_cast(&bounds[0]);
	// Instantiate the random number adjuster
	adjust_randoms_functor adjuster(cumul_prob_bounds);

	// Adjust the random numbers to fall inside the correct interval
	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(
				uniform_begin, 
				inds_demes_begin
				)),
			thrust::make_zip_iterator(thrust::make_tuple(
				uniform_end,
				inds_demes_end					
				)),
				adjuster
				);
	}