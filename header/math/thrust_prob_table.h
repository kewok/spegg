#ifndef THRUST_PROB_TABLE_H
#define THRUST_PROB_TABLE_H

#include <thrust/host_vector.h>

class ThrustProbTable
{
/* 
* A look-up table to be used for simulating from arbitrary discrete distributions.
*/
	public:
		void setup(thrust::host_vector<float>::iterator prob_begin, thrust::host_vector<float>::iterator prob_end);
		void draw(thrust::host_vector<float>::iterator uniform_begin, thrust::host_vector<float>::iterator uniform_end, thrust::host_vector<int>::iterator result);
	protected:
		thrust::host_vector<float> cumulative_prob;
};

#endif
