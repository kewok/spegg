#ifndef MATING_THRUST_PROB_TABLE_DEMES_H
#define MATING_THRUST_PROB_TABLE_DEMES_H

#include <math/thrust_prob_table.h>

class mating_ThrustProbTable_demes : public ThrustProbTable
	{
	public:
		void adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);

		void adjust_randoms_fixed_offsets(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);

		void determine_key_offsets(int number_of_key_types, thrust::device_vector<int> &key_histogram_vector);

		thrust::device_vector<int> key_offsets;
	};

class mating_ThrustProbTable_demes_Double : public ThrustProbTableDouble
	{
	public:
		void adjust_randoms(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);

		void adjust_randoms_fixed_offsets(thrust::device_vector<double>::iterator uniform_begin, thrust::device_vector<double>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);

		void determine_key_offsets(int number_of_key_types, thrust::device_vector<int> &key_histogram_vector);

		thrust::device_vector<int> key_offsets;
	};

#endif
