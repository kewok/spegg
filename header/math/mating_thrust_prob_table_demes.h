#ifndef MATING_THRUST_PROB_TABLE_DEMES_H
#define MATING_THRUST_PROB_TABLE_DEMES_H

#include <math/thrust_prob_table.h>

class mating_ThrustProbTable_demes : public ThrustProbTable
	{
	public:
		void adjust_randoms(thrust::host_vector<float>::iterator uniform_begin, thrust::host_vector<float>::iterator uniform_end,
		thrust::host_vector<int>::iterator inds_demes_begin, thrust::host_vector<int>::iterator inds_demes_end);

		void adjust_randoms_fixed_offsets(thrust::host_vector<float>::iterator uniform_begin, thrust::host_vector<float>::iterator uniform_end, thrust::host_vector<int>::iterator inds_demes_begin, thrust::host_vector<int>::iterator inds_demes_end);

		void determine_key_offsets(int number_of_key_types, thrust::host_vector<int>& key_histogram_vector );

		thrust::host_vector<int> key_offsets;
	};

#endif
