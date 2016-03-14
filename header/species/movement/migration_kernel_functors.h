#ifndef MIGRATION_KERNEL_FUNCTORS_H
#define MIGRATION_KERNEL_FUNCTORS_H

#include <species/inds.h>
#include "math.h"
#include <util/thrust_functors.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>

/* Functors related to migration */

struct determine_if_migratory
	{
/* 
* A thrust functor that draws a Bernoulli random variable which determines whether an individual migrates during each time step. The migration probability is presumably calculated separately, as is the uniform random number. 
*/
	/* 
		Elements in the tuple.


		----------------------
		0: whether the individual migrates
		1: their migration probability
		2: the associated random number 
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<2>(t) < thrust::get<1>(t))
			thrust::get<0>(t) = 1;
		}
	};


struct move_individuals_functor
	{
/* 
* A thrust functor that simulates migration by reassigning the deme of the migrating individual. 
*/

	/* 
		Elements in the tuple.


		----------------------
		0: individual's deme
		1: whether the individual migrates
		2: individual's destination
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<0>(t) = thrust::get<1>(t)*thrust::get<2>(t) + (1 - thrust::get<1>(t))*thrust::get<0>(t);
		}
	};
#endif
