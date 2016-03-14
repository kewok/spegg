#ifndef SAMPLING_EVENT_H
#define SAMPLING_EVENT_H

#include <iostream>
#include <stdio.h>
#include <thrust/host_vector.h>

#include <util/sampling_input.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class SamplingEvent
	{
	public:
	/* use a factory method */
		static SamplingEvent *create_SamplingEvent(SamplingInput *sampling_input, gsl_rng* gen);
	
		thrust::host_vector<int> sampling_individuals_indices;
		thrust::host_vector<int> sampled_individuals_indices;

	virtual void sample() = 0;
	
	protected:
		SamplingInput *sampling_input;
		gsl_rng* gen;
	};

#endif
