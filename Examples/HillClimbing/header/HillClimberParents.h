#ifndef HILLCLIMBER_PARENTS_CLASS_H
#define HILLCLIMBER_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>

class HillClimberParents : public Parents
	{
	class HillClimbers *species;
	public:
		HillClimberParents(HillClimbers *species);
		
	protected:
		int FECUNDITY_PHENOTYPE_INDEX;

		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();
	};


struct reproductive_probability_functor
	{
	float *reproductive_advantage;

	reproductive_probability_functor(float* reproductive_adv) : reproductive_advantage(reproductive_adv)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: whether the individual will reproduce
		1: the individual's deme
		2: the individual's fecundity value
		3: the probability the individual becomes a parent.
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<3>(t) = thrust::get<0>(t) * powf(thrust::get<2>(t), reproductive_advantage[thrust::get<1>(t)]);
		}
	};

#endif 
