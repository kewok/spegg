#include <curand.h>
#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
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
#include <thrust/extrema.h>

#include <util/thrust_functors.h>
#include <math/thrust_probabilities.h>
#include <math/histogram.h>
/* 

Stolen from:
http://code.google.com/p/thrust/source/browse/examples/histogram.cu

*/

// TODO: Add a case for handling what to do when data are of size zero.

void calculate_histogram(thrust::device_vector<int> &data,
			 thrust::device_vector<int> &counts,
			 int counts_size)
	{
	counts.resize( counts_size );
	
	int data_size = data.size();
	thrust::counting_iterator<int> search_begin(0);
	thrust::device_vector<int> temp_data( data_size );

	thrust::copy(data.begin(), data.end(), temp_data.begin());

	thrust::sort(temp_data.begin(), temp_data.end());

	thrust::upper_bound(temp_data.begin(),  temp_data.end(), 
                      search_begin, search_begin + counts_size,
                      counts.begin());

	thrust::adjacent_difference(counts.begin(), counts.end(),
                              counts.begin());
	}


void calculate_histogram(thrust::device_vector<float> &data,
			 thrust::device_vector<int> &counts,
			 int counts_size)
	{
	counts.resize( counts_size );
	
	int data_size = data.size();
	thrust::counting_iterator<int> search_begin(0);
	thrust::device_vector<float> temp_data( data_size );

	float maxval=*(thrust::max_element(data.begin(), data.end()));
	float minval=*(thrust::min_element(data.begin(), data.end()));

	thrust::device_vector<float> maxes( data_size );
	thrust::device_vector<float> mins( data_size );
	thrust::device_vector<int> bin_nums( data_size );

	thrust::fill(maxes.begin(), maxes.end(), maxval);
	thrust::fill(mins.begin(), mins.end(), minval);
	thrust::fill(bin_nums.begin(), bin_nums.end(), counts_size);
	
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(data.begin(), mins.begin(), maxes.begin(), bin_nums.begin(), temp_data.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(data.end(), mins.end(), maxes.end(), bin_nums.end(), temp_data.end())),
			 adjust_histogram_data_values());

	thrust::sort(temp_data.begin(), temp_data.end());

	thrust::upper_bound(temp_data.begin(),  temp_data.end(), 
                      search_begin, search_begin + counts_size,
                      counts.begin());

	thrust::adjacent_difference(counts.begin(), counts.end(),
                              counts.begin());	
	}


/*****************
// To extract the counts for a subset of the histogram elements pertaining to specific values
// e.g., data=(0,0,0,1,1,2,2,2,2,3) with histogram counts = (3,2,4,1), extract only counts for element 1 and 3 so end up with histogram = (2,1)

// The variable counts_size refers to the number of values data can take on

Todo: Add an error catch that if counts_size < max(data), return error 
//
*****************/

void calculate_histogram_subset(thrust::device_vector<int> &data,
			 thrust::device_vector<int> &counts,
			 int counts_size,
			 thrust::device_vector<int> &subset_indices)
	{
	calculate_histogram(data, counts, counts_size);

	thrust::device_vector<int> final_counts( subset_indices.size() );

	thrust::gather(subset_indices.begin(), subset_indices.end(), counts.begin(), final_counts.begin());

	counts.resize( subset_indices.size() );

	thrust::copy( final_counts.begin(), final_counts.end(), counts.begin() );
	}

