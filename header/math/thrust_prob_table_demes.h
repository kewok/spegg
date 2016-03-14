#ifndef THRUST_PROB_TABLE_DEMES_H
#define THRUST_PROB_TABLE_DEMES_H

#include <math/thrust_prob_table.h>

class ThrustProbTable_demes : public ThrustProbTable
	{
	/* 
	* A derived version of the ThrustProbTable capable of adjusting the look-up probabilities to fall within ranges specified by the deme: e.g., if you have individuals with demes [0,0,1,1] and look up probabilities [0.5,0.5, 0.5, 0.5], a ThrustProbTable_demes object transforms this into [0.25,0.5,0.75,1.0], and the random draws in deme 0 will be from U(0,0.5) and for deme 1 from U(0.5,1) so that deme-specific random draws can be made in parallel.
	*/

	public:
		void adjust_randoms(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end,
		thrust::device_vector<int>::iterator deme_offsets_begin, thrust::device_vector<int>::iterator deme_offsets_end,
		thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);

		void adjust_randoms_fixed_offsets(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator inds_deme_begin, thrust::device_vector<int>::iterator inds_deme_end);
	};

#endif
