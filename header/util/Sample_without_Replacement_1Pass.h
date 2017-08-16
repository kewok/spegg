#ifndef SAMPLE_WITHOUT_REPLACEMENT_1PASS_H
#define SAMPLE_WITHOUT_REPLACEMENT_1PASS_H

#include <iostream>
#include <stdio.h>
#include <util/sampling_event.h>

class Sample_without_Replacement_1Pass : public SamplingEvent
	{
	public:
		Sample_without_Replacement_1Pass(SamplingInput *sampling_input, curandGenerator_t gen )
			{
			this->sampling_input = sampling_input;
			this->gen = gen;
			}
	
	void sample();
	};

#endif
