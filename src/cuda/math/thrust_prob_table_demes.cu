#include <math/thrust_prob_table_demes.h>
#include <util/footimer2.h>
#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/gather.h>

// make sure the individual's uniform random number falls in the interval corresponding to their current deme
struct adjust_randoms_functor
	{
	float *popRanges_in_cumulative;

	adjust_randoms_functor(float *pop_range) : popRanges_in_cumulative(pop_range)
	{};
	
	/*
		Elements in the tuple.
		----------------------
		0: the individual's random number
		1: the individual's deme
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
			int this_inds_deme = thrust::get<1>(t);
			if (this_inds_deme != 0)
				{
				float interval = popRanges_in_cumulative[this_inds_deme] - popRanges_in_cumulative[this_inds_deme - 1];
				thrust::get<0>(t) = thrust::get<0>(t)*interval +  popRanges_in_cumulative[this_inds_deme - 1];
				}

	 // For the case of the individuals in deme zero
			else
				{
				float interval = popRanges_in_cumulative[this_inds_deme];
				thrust::get<0>(t) = thrust::get<0>(t)*interval;
				}
		}
	};

struct adjust_randoms_functor_double
	{
	double *popRanges_in_cumulative;

	adjust_randoms_functor_double(double *pop_range) : popRanges_in_cumulative(pop_range)
	{};
	
	/*
		Elements in the tuple.
		----------------------
		0: the individual's random number
		1: the individual's deme
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
			int this_inds_deme = thrust::get<1>(t);
			if (this_inds_deme != 0)
				{
				double interval = popRanges_in_cumulative[this_inds_deme] - popRanges_in_cumulative[this_inds_deme - 1];
				thrust::get<0>(t) = thrust::get<0>(t)*interval +  popRanges_in_cumulative[this_inds_deme - 1];
				}

	 // For the case of the individuals in deme zero
			else
				{
				double interval = popRanges_in_cumulative[this_inds_deme];
				thrust::get<0>(t) = thrust::get<0>(t)*interval;
				}
		}
	};


void ThrustProbTable_demes::adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end,
thrust::device_vector<int>::iterator deme_offsets_begin, thrust::device_vector<int>::iterator deme_offsets_end,
thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end)
	{
	// Determine the number of demes.
	int n = thrust::distance(deme_offsets_begin, deme_offsets_end);

	// Figure out the bounds for each deme in terms of their values in the cumulative probability table
	thrust::device_vector<float> bounds(n);

	thrust::gather(deme_offsets_begin, deme_offsets_end,cumulative_prob.begin(),bounds.begin());
	
	thrust::device_vector<int> temp(n);
	thrust::copy(deme_offsets_begin, deme_offsets_end, temp.begin());

	float *cumul_prob_bounds = raw_pointer_cast(&bounds[0]);
	// Instantiate the random number adjuster
	adjust_randoms_functor adjuster(cumul_prob_bounds);

	// Adjust the random numbers to fall inside the correct interval
	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(
				uniform_begin, 
				inds_deme_begin
				)),
			thrust::make_zip_iterator(thrust::make_tuple(
				uniform_end,
				inds_deme_end 
				)),
				adjuster
				);
	}



void ThrustProbTable_demes_Double::adjust_randoms(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end,
thrust::device_vector<int>::iterator deme_offsets_begin, thrust::device_vector<int>::iterator deme_offsets_end,
thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end)
	{
	// Determine the number of demes.
	int n = thrust::distance(deme_offsets_begin, deme_offsets_end);

	// Figure out the bounds for each deme in terms of their values in the cumulative probability table
	thrust::device_vector<double> bounds(n);

	thrust::gather(deme_offsets_begin, deme_offsets_end,cumulative_prob.begin(),bounds.begin());
	
	thrust::device_vector<int> temp(n);
	thrust::copy(deme_offsets_begin, deme_offsets_end, temp.begin());

	double *cumul_prob_bounds = raw_pointer_cast(&bounds[0]);
	// Instantiate the random number adjuster
	adjust_randoms_functor_double adjuster(cumul_prob_bounds);

	// Adjust the random numbers to fall inside the correct interval
	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(
				uniform_begin, 
				inds_deme_begin
				)),
			thrust::make_zip_iterator(thrust::make_tuple(
				uniform_end,
				inds_deme_end 
				)),
				adjuster
				);
	}

