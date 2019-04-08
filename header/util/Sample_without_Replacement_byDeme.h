#ifndef SAMPLE_WITHOUT_REPLACEMENT_ByDeme_H
#define SAMPLE_WITHOUT_REPLACEMENT_ByDeme_H

#include <thrust/sort.h>
#include <thrust/copy.h>
#include <util/sampling_event.h>

// This version can be used to stratify sampling by deme, provided the individuals potentially conducting the sampling only sample one individual in each time step.

class Sample_without_Replacement_byDeme : public SamplingEvent
	{
	public:
		Sample_without_Replacement_byDeme(SamplingInput *sampling_input, curandGenerator_t gen)
			{
			this->sampling_input = sampling_input;
			this->gen = gen;

			sampling_individuals_indices.resize(sampling_input->list_of_individuals_potentially_conducting_sampling.size());
			thrust::copy(sampling_input->list_of_individuals_potentially_conducting_sampling.begin(), sampling_input->list_of_individuals_potentially_conducting_sampling.end(), sampling_individuals_indices.begin());

			sampled_individuals_indices.resize(sampling_input->list_of_individuals_potentially_conducting_sampling.size());

			index_to_sample.resize(sampling_individuals_indices.size());
			unique_uniform_rvs.resize(sampling_input->list_of_individuals_potentially_subject_to_sampling.size());
			demes_of_individuals_subject_to_sampling.resize(sampling_input->list_of_individuals_potentially_subject_to_sampling.size());

			cumulative_sampleable_individuals_by_deme.resize(sampling_input->Num_Demes);
			number_of_sampling_individuals_by_deme.resize(sampling_input->Num_Demes);
			cumulative_sampling_individuals_by_deme.resize(sampling_input->Num_Demes);
			}
		void sample();
		void setup_demes(thrust::device_vector<int> &demes_of_individuals_sampled);

	protected:
		thrust::device_vector<int> cumulative_sampleable_individuals_by_deme;
		thrust::device_vector<int> cumulative_sampling_individuals_by_deme;
		thrust::device_vector<int> number_of_sampling_individuals_by_deme;
		thrust::device_vector<double> unique_uniform_rvs;
		thrust::device_vector<int> demes_of_individuals_subject_to_sampling;
		thrust::device_vector<int> index_to_sample;
	};


// A functor which, for each indiviudal conducting the sampling, assigns the index of the individuals subject to sampling which the individual conducting the sampling will pick. This is only used when sampling is stratified by demes. Note this assumes that the number of individuals conducting the sampling is less than or equal to the number of individuals being sampled.

struct specify_index_to_sample
	{
	int *cumulative_demewise_sum_individuals_conducting_sampling;
	int *cumulative_demewise_sum_individuals_subject_to_sampling;

	specify_index_to_sample(int *cumsum_demes_focal, int *cumsum_demes_target) : cumulative_demewise_sum_individuals_conducting_sampling(cumsum_demes_focal), cumulative_demewise_sum_individuals_subject_to_sampling(cumsum_demes_target)
	{};
	/*
		Elements in the tuple.
		----------------------
		0: index of individuals conducting the sampling
		1: deme of the individuals conducting the sampling
		2: the indices of the individuals subject to sampling by the individuals conducting the sampling (should be of same length as indices of individuals conducting the sampling)
	*/

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ind_deme = thrust::get<1>(t);
		int ind_val = thrust::get<0>(t);

		thrust::get<2>(t) = thrust::get<0>(t);
		int ans = 0;
		if (ind_deme > 0)
			{
			ans = thrust::get<0>(t)% (cumulative_demewise_sum_individuals_conducting_sampling[ind_deme]);
			thrust::get<2>(t) = ans + (cumulative_demewise_sum_individuals_subject_to_sampling[ind_deme]); 
			}
		}
	};

#endif
