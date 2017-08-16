#ifndef SAMPLING_EVENT_H
#define SAMPLING_EVENT_H

#include <iostream>
#include <stdio.h>
#include <util/sampling_input.h>


class SamplingEvent
	{
	public:
	/* use a factory method */
		static SamplingEvent *create_SamplingEvent(SamplingInput *sampling_input, curandGenerator_t gen);
	
		thrust::device_vector<int> sampling_individuals_indices;
		thrust::device_vector<int> sampled_individuals_indices;

		virtual void sample() = 0;
	
	protected:
		SamplingInput *sampling_input;
		curandGenerator_t gen;
	};

#endif
