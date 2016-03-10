#ifndef THRUST_PROBABILITIES_H
#define THRUST_PROBABILITIES_H

#define PI 3.14159265358979f
#define POISSON_MAX_RV 0.99999991f
#define MAXIMUM_ITERATIONS_REJECTION_SAMPLING 1000000

#include "math.h"
#include <thrust/functional.h>
#include <stdio.h>

// Takes the float value given by curandGenerateUniform() and converts it into a discrete value between 0, n-1
struct discrete_uniform
	{
	const int maxSize;
	discrete_uniform(int _maxSize) : maxSize(_maxSize) {};

	__host__ __device__
	int operator()(const float& rv) const {
		int ans = 0;
		float ansf = ((float) maxSize)*rv;
		ans = (int) ansf;
		return ans;
		}
	};

struct bernoulli_rv
	{
	const float probability;
	bernoulli_rv(float _probability) : probability(_probability) {};
	
	__host__ __device__
	int operator()(const float& rv) const {
		int ans = 0;
		if (rv <= probability)
			ans = 1;
		if (rv > probability)
			ans = 0;
		return ans;
		}
	};


struct poisson_rv
	{
	const float lambda;
	poisson_rv(float _lambda) : lambda(_lambda) {};
	
	__host__ __device__
	int operator()(const float& rv) const {
		float use_rv = rv;
		if (use_rv > POISSON_MAX_RV)
			{
			use_rv = POISSON_MAX_RV;
			}
		int ans = 0;
		float p = exp(-lambda);
		float F = p;
		int maxits = 0;
		while ((use_rv > F) && (maxits < MAXIMUM_ITERATIONS_REJECTION_SAMPLING))
			{
			p = lambda*p/(ans + 1);
			F = F + p;
			ans = ans + 1;
			maxits++;
			}
		return ans;
		}
	};

struct discrete_normal_rv
	{
	const float mean;
	const float sd;

	discrete_normal_rv(float _mean, float _sd) : mean(_mean), sd(_sd) {};

	/*
		Elements in the tuple.
		----------------------
		0: uniform rv 1
		1: uniform rv 2
		2: ans
	*/

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<2>(t) = (int) sd * r * cosf(phi) + mean;
		}
	};

struct normal_rv
	{
	const float mean;
	const float sd;

	normal_rv(float _mean, float _sd) : mean(_mean), sd(_sd) {};

	/*
		Elements in the tuple.
		----------------------
		0: uniform rv 1
		1: uniform rv 2
		2: ans
	*/

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<2>(t) = sd * r * cosf(phi) + mean;
		}
	};

struct normal_rv_different_parameters
	{

	/*
		Elements in the tuple.
		----------------------
		0: uniform rv 1
		1: uniform rv 2
		2: mean
		3: sd
		4: ans
	*/

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<4>(t) = thrust::get<3>(t) * r * cosf(phi) + thrust::get<2>(t);
		}
	};


struct discrete_normal_rv_different_parameters
	{

	/*
		Elements in the tuple.
		----------------------
		0: uniform rv 1
		1: uniform rv 2
		2: mean
		3: sd
		4: ans
	*/

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float r = sqrtf(-2.0f * logf(thrust::get<0>(t)));
		float phi = 2 * PI * thrust::get<1>(t);
		thrust::get<4>(t) = (int) (thrust::get<3>(t) * r * cosf(phi) + thrust::get<2>(t));
		}
	};

struct poisson_rv_different_parameters
	{
	/* 
		Elements in the tuple.
		----------------------
		0: uniform rv
		1: lambda
		2: ans
	*/	

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ans = 0;
		float lambda = thrust::get<1>(t);
		float rv = thrust::get<0>(t);
		if (rv > POISSON_MAX_RV)
			{
			rv = POISSON_MAX_RV;
			}

		float p = exp(-lambda);
		float F = p;
		int maxits = 0;
		while ((rv > F) && (maxits < MAXIMUM_ITERATIONS_REJECTION_SAMPLING))
			{
			maxits++;
			p = lambda*p/(ans + 1);
			F = F + p;
			ans = ans + 1;
			}
		thrust::get<2>(t) = ans;
		}
	};

struct bernoulli_rv_different_parameters
	{
	/* 
		Elements in the tuple.
		----------------------
		0: uniform rv
		1: probability of success
		2: ans
	*/	

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ans = 0;
		float prob = thrust::get<1>(t);
		float rv = thrust::get<0>(t);
		if (rv <= prob)
			{
			ans = 1;
			}
		if (rv > prob)
			{
			ans = 0; 
			}
		thrust::get<2>(t) = ans;
		}
	};
#endif
