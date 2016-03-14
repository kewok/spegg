#ifndef PREY_H
#define PREY_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <util/footimer2.h>
#include <util/thrust_functors.h>
#include "prey_deme_specific_data.h"

class prey
	{
	public:
		prey(int num_demes, int prey_index, int seed_val, int num_time_steps);
		~prey();

		// Random number generator
		curandGenerator_t gen;
	
		// Misc data integers
		int seed;
		int Num_Demes;
		int Prey_Index;
		int Number_of_Timesteps;

		// Data vectors, each of length Num_Demes, keeping track of different data relevant to the prey dynamics
		thrust::device_vector<float> prey_abundance;
		thrust::device_vector<float> effect_of_inds_on_prey;
		thrust::device_vector<float> prey_maximum_abundance;
		thrust::device_vector<float> prey_assimilation_efficiency;

		void specify_prey_properties_by_deme();

		void update_prey_abundance(thrust::device_vector<float> &effect_of_individuals_on_prey);

	protected:
		PreyDemeSpecificData *preyParameters;

		thrust::device_vector<float> prey_carrying_capacity;
		thrust::device_vector<float> prey_unconstrained_growth_rates;
		thrust::device_vector<float> prey_growth_rate_noise_stddev;

		// Functions related to initializations

		void calculate_maximum_abundance();
	};


struct update_prey
{
	/* 

		Elements in the tuple.

		----------------------
		0: prey unconstrained growth rate
		1: prey carrying capacity
		2: effect of inds on prey
		3. prey abundance
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		thrust::get<3>(t) = thrust::get<0>(t)*thrust::get<3>(t)/(1+thrust::get<3>(t)/thrust::get<1>(t)) - thrust::get<2>(t); 
		if (thrust::get<3>(t) < 0)
			{
			thrust::get<3>(t) = 0; // Lower ceiling
			}
		}	
};

#endif 
