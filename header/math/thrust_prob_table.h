#ifndef THRUST_PROB_TABLE_H
#define THRUST_PROB_TABLE_H

#include <thrust/device_vector.h>

class ThrustProbTable
	{
/* 
* A look-up table to be used for simulating from arbitrary discrete distributions.
*/
	public:
		void setup(thrust::device_vector<float>::iterator prob_begin, thrust::device_vector<float>::iterator prob_end);
		void draw(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator result);
		thrust::device_ptr<float> get_cumulative_probabilities_ptr();

	protected:
		thrust::device_vector<float> cumulative_prob;
	};

#endif
