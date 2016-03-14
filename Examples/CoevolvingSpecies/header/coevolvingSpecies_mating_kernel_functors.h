#ifndef COEVOLVING_SPECIES_MATING_KERNEL_FUNCTORS_H
#define COEVOLVING_SPECIES_MATING_KERNEL_FUNCTORS_H

struct coevolving_parental_eligibility_functor
{
	int parental_sex;
	float *age_at_maturity;

	coevolving_parental_eligibility_functor(float* age_mat, int sex) : parental_sex(sex), age_at_maturity(age_mat)
	{};

	/* 
		Elements in the tuple.

		----------------------

		0: whether the individual will reproduce
		1: the individual's sex
		2: the individual's deme
		3: the individual's age


	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<1>(t)==parental_sex)
			{
			float ind_age = (float) thrust::get<3>(t);
			if (age_at_maturity[thrust::get<2>(t)] <= ind_age)
				thrust::get<0>(t) = 1;
			}
		}
};

struct bisexual_parental_eligibility_functor
{
	float *age_at_maturity;

	bisexual_parental_eligibility_functor(float* age_mat) : age_at_maturity(age_mat)
	{};

	/* 
		Elements in the tuple.

		----------------------

		0: whether the individual will reproduce
		1: the individual's deme
		2: the individual's age


	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float ind_age = (float) thrust::get<2>(t);
			
		if (age_at_maturity[thrust::get<1>(t)]<=ind_age)
			thrust::get<0>(t) = 1;
		}
};

#endif

