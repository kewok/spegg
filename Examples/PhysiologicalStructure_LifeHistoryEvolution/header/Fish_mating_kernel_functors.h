#ifndef FISH_MATING_KERNEL_FUNCTORS_H
#define FISH_MATING_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/host_vector.h>
#include <thrust/functional.h>

struct parental_eligibility_functor
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
	__host__
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
	__host__ 
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
	__host__
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

#endif
