#ifndef THRUST_FUNCTORS_H
#define THRUST_FUNCTORS_H

#include <thrust/tuple.h>

#define MAX_RV 0.999991f // Used because floats occasionally cause indexing problems in reassign_functor().

template <typename T>
struct unary_plus
	{
	T rhs;
	unary_plus(const T &r) : rhs(r) {};
	
	__host__ __device__
	T operator()(const T &lhs) {
		return lhs + rhs;
	}
	};

template <typename T>
struct unary_minus
	{
	T rhs;
	unary_minus(const T &r) : rhs(r) {};
	
	__host__ __device__
	T operator()(const T &lhs) {
		return lhs - rhs;
	}
	};

template <typename T>
struct unary_modulus
	{
	T rhs;
	unary_modulus(const T &r) : rhs(r) {};
	
	__host__ __device__
	T operator()(const T &lhs) {
		return lhs%rhs;
	}
	};

template <typename T>
struct unary_less
	{
	T rhs;
	unary_less(const T &r) : rhs(r) {};
	
	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs < rhs;
	}
	};

template <typename T>
struct unary_less_equal
	{
	T rhs;
	unary_less_equal(const T &r) : rhs(r) {};
	
	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs <= rhs;
	}
	};

template <typename T>
struct unary_greater
	{
	T rhs;
	unary_greater(const T &r) : rhs(r) {};
	
	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs > rhs;
	}
	};

template <typename T>
struct unary_not_equal
	{
	T rhs;
	unary_not_equal(const T &r) : rhs(r) {};
	
	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs != rhs;
	}
	};

template <typename T>
struct unary_greater_equal
	{
	T rhs;
	unary_greater_equal(const T &r) : rhs(r) {};
	
	__host__ __device__
	bool operator()(const T &lhs) {
		return lhs >= rhs;
	}
	};

struct is_less_than_zero
  {
    __host__ __device__
    bool operator()(int x)
    {
      return x < 0;
    }
  };


struct is_less_than_zero_f
  {
    __host__ __device__
    bool operator()(float x)
    {
      return x < 0;
    }
  };

struct pairwise_min
	{
	/*
	Elements in the tuple.
	----------------------
	0: LHS
	1: RHS
	2: Pairwise minimum value
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
	if (thrust::get<0>(t) <= thrust::get<1>(t))
		thrust::get<2>(t) = thrust::get<0>(t);
	else
		thrust::get<2>(t) = thrust::get<1>(t);
	}
	};

/* The functor reassign_functor takes a float and scales it to fall within a given integer interval */

/* The functor reassign_functor takes a number and scales it to fall within a given interval */

struct reassign_functor
	{
	int *sampled_offset;
	reassign_functor(int *d_offset) : sampled_offset(d_offset)
	{};

        /*
	Elements in the tuple.
	----------------------
	0: the sampled individuals's floating point uniform random number
	1: the population from which the individual doing the sampling comes
	2: the sampled individual's value
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
	int sampler_from = thrust::get<1>(t);
	
 // Assign the sampled individual to the new id accordingly
	if (sampler_from != 0)	
		{
		thrust::get<2>(t) = (int) (thrust::get<0>(t)*(sampled_offset[sampler_from] - sampled_offset[sampler_from-1]) + (sampled_offset[sampler_from-1]));

		// Apply correction caused by floating point approximation issues. 

		if (thrust::get<2>(t) == (sampled_offset[sampler_from] - sampled_offset[sampler_from-1]) + (sampled_offset[sampler_from-1]) )
			{
			thrust::get<2>(t) = thrust::get<2>(t)-1; 
			}
		}
	else
		{
		thrust::get<2>(t) = (int) (thrust::get<0>(t)*(sampled_offset[sampler_from]));
		if (thrust::get<2>(t) == sampled_offset[sampler_from])
			{
			thrust::get<2>(t) = sampled_offset[sampler_from]-1; 
			}
	
		}
	}
	};

struct reassign_functor_simplified
	{
        /*
	Elements in the tuple.
	----------------------
	0: the floating point number (e.g., r.v~U(0,1))
	1: the lower value of the integer interval -> determining these values can be quite tricky and inefficient
	2: the upper value of the integer interval
	3: the integer value answer
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
	if (thrust::get<0>(t) > MAX_RV)
			{
			thrust::get<0>(t) = MAX_RV;
			}

	thrust::get<3>(t) = (int) (thrust::get<0>(t) * ( thrust::get<2>(t) - thrust::get<1>(t) ));
	}
	};


// functor that changes floats to ints
struct float_to_int
	{
	/*
		Elements in the tuple.
		----------------------
		0: float input
		1: int output
	*/
	template <typename tuple>
	__host__ __device__ 
	void operator() ( tuple t ) {
		thrust::get<1>(t) = floorf(thrust::get<0>(t) + 0.5);
		}
	};

struct int_to_float
	{
/*
* custom functor that changes floats to ints
/
	/*
		Elements in the tuple.
		----------------------
		0: int input
		1: float output
	*/
	template <typename tuple>
	__host__ __device__ 
	void operator() ( tuple t ) {
		thrust::get<1>(t) = (float) thrust::get<0>(t);
		}
	};

struct int_to_double
	{
/*
* custom functor that changes floats to ints
/
	/*
		Elements in the tuple.
		----------------------
		0: int input
		1: double output
	*/
	template <typename tuple>
	__host__ __device__ 
	void operator() ( tuple t ) {
		thrust::get<1>(t) = (double) thrust::get<0>(t);
		}
	};

// functor that returns absolute value of a float
struct fabs_functor
	{
	/*
		Elements in the tuple.
		----------------------
		0: float value
	*/
	template <typename tuple>
	__host__ __device__ 
	void operator() ( tuple t ) {
		if (thrust::get<0>(t) < 0)
			{
			thrust::get<0>(t) = -1*thrust::get<0>(t);
			}
		}
	};


#endif
