#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/replace.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/unique.h>
#include <thrust/set_operations.h>

#include <util/thrust_functors.h>
#include <math/thrust_probabilities.h>


void calculate_histogram(thrust::host_vector<int> &data, thrust::host_vector<int> &counts, 		 int counts_size);

void calculate_histogram(thrust::host_vector<float> &data, thrust::host_vector<int> &counts, 		 int counts_size);

void calculate_histogram_subset(thrust::host_vector<int> &data, thrust::host_vector<int> &counts, int counts_size, thrust::host_vector<int> &subset_indices);


// convert values to fall inside bin groups
struct adjust_histogram_data_values
{
	/* 
		Elements in the tuple.
		----------------------
		0: original data value
		1: minimum data value
		2: maximum data value
		3: number of bins
		4: new data value
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<4>(t) = thrust::get<0>(t) - thrust::get<1>(t);
		thrust::get<4>(t) = thrust::get<4>(t)/(thrust::get<2>(t) - thrust::get<1>(t));
		thrust::get<4>(t) *= (float) thrust::get<3>(t); 
		}
	
};

#endif
