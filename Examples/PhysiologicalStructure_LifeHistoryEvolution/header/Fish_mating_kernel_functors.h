#ifndef FISH_MATING_KERNEL_FUNCTORS_H
#define FISH_MATING_KERNEL_FUNCTORS_H

#include <thrust/functional.h>

typedef thrust::device_vector<float>::iterator floatIter_t;
typedef thrust::device_vector<int>::iterator intIter_t;

struct parental_eligibility_functor
/* note this overwrites the parental_eligibility_functor in mating_kernel_functors.h */
{
	int parental_sex;
	float *size_at_maturity;

	parental_eligibility_functor(float* size_mat, int sex) : parental_sex(sex), size_at_maturity(size_mat)
	{};

	/* 
		Elements in the tuple.

		----------------------

		0: whether the individual will reproduce
		1: the individual's sex
		2: the individual's subpopulation
		3: the individual's size


	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<1>(t)==parental_sex)
			{
			if (size_at_maturity[thrust::get<2>(t)]<=thrust::get<3>(t))
				thrust::get<0>(t) = 1;
			}
		}
};

struct reproductive_probability_functor
{
	int parental_sex;
	float *reproductive_advantage;

	reproductive_probability_functor(float* reproductive_adv, int sex) : parental_sex(sex), reproductive_advantage(reproductive_adv)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: whether the individual will reproduce
		1: the individual's subpopulation
		2: the number of kids the individual will have
		3: the individual's reproductive probability

	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float num_kids_by_mom = (float) thrust::get<2>(t);
		thrust::get<3>(t) = thrust::get<0>(t) * powf(num_kids_by_mom, reproductive_advantage[thrust::get<1>(t)]);
		}
};

struct female_fecundity_functor
{
	/* 
		Elements in the tuple.

		----------------------
		0: whether the individual is a reproductive female
		1: the individual's fecundity phenotype
		2: the individual's final fecundity score
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<0>(t) > 0)
			{
			thrust::get<2>(t) = thrust::get<1>(t) ;
			}

		// This hack is added to prevent integer overflow. In a later version of this code, int should be replaced with unsigned int
		if (thrust::get<2>(t) > 5000)
			{
			thrust::get<2>(t) = 5000;
			}
		}
};

struct determine_female_nesting_success_functor
{
	float *number_of_available_nests;
	float *number_of_reproductive_females; 
	int *number_of_reproductive_males;
	
	determine_female_nesting_success_functor(float *number_of_available_nests_ptr, float *number_of_reproductive_females_ptr, int* number_of_reproductive_males_ptr) : number_of_available_nests(number_of_available_nests_ptr), number_of_reproductive_females(number_of_reproductive_females_ptr), number_of_reproductive_males(number_of_reproductive_males_ptr)
	{};

	/* 
		Elements in the tuple.
pop_begin, will_reproduce_begin, probability_female_nests.begin(), rand.begin()

		----------------------
		0: whether the individual will reproduce
		1: the individual's subpopulation
		2: the individual's random number
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<0>(t) == 1) /* if the individual is indeed a reproducing female */
			{
			int individuals_subpopulation = thrust::get<1>(t);
		
			float nesting_probability = number_of_available_nests[individuals_subpopulation]/(1+number_of_reproductive_females[individuals_subpopulation]); // The 1+x is needed to prevent cases when there are no reproductive females

			if (thrust::get<2>(t) <= nesting_probability)
				{
				thrust::get<0>(t) = 1;
				}
			else
				{
				thrust::get<0>(t) = 0;
				}

			if (number_of_reproductive_males[individuals_subpopulation] == 0)
				{
				thrust::get<0>(t) = 0;
				}

			}
		}
};
#endif
