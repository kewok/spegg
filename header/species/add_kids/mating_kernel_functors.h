#ifndef REPRODUCTION_KERNEL_FUNCTORS_H
#define REPRODUCTION_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>


struct parental_eligibility_functor
	{
	int parental_sex;

	parental_eligibility_functor(int sex) : parental_sex(sex)
	{};

	/* 
		Elements in the tuple.

		----------------------

		0: whether the individual will reproduce
		1: the individual's sex

	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<1>(t)==parental_sex)
			{
			thrust::get<0>(t) = 1;
			}
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

		}
	};

#endif
