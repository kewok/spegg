#ifndef SAMPLE_WITHOUT_REPLACEMENT_Test_H
#define SAMPLE_WITHOUT_REPLACEMENT_Test_H

#include <thrust/sort.h>
#include <thrust/copy.h>
#include <util/sampling_event.h>

class Sample_without_Replacement_Test : public SamplingEvent
	{
	public:
		Sample_without_Replacement_Test(SamplingInput *sampling_input, curandGenerator_t gen )
			{
			this->sampling_input = sampling_input;
			this->gen = gen;

			number_of_individuals_subject_to_sampling = sampling_input->list_of_individuals_potentially_subject_to_sampling.size();
			}
		void sample();

	protected:
		int number_of_individuals_subject_to_sampling;
		thrust::device_vector<double> unique_uniform_rvs;
		void draw_unique_double_randoms();
		void shuffle_sampled_individuals_at_random();	
	};

#endif
