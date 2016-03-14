#ifndef INTERACTION_KERNEL_FUNCTORS_H
#define INTERACTION_KERNEL_FUNCTORS_H

#include <curand.h>
#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/random.h>

#define PI 3.14159265358979f
#define POISSON_MAX_RV 0.99999991f
#define MAXIMUM_ITERATIONS_REJECTION_SAMPLING 1000000
#define MAX_SAMPLED 20

inline __host__ __device__ int draw_ind(int min_ind_index, int max_ind_index, float use_rv)
{
	
	int range = max_ind_index - min_ind_index;
	int ans = ((int) (use_rv*((float) range))) + min_ind_index;
	return ans;
}

inline __host__ __device__ float calculate_effect_of_interaction(float effect_of_interaction, float weighted_effect_of_alternative_individual, float interaction_strength, float focal_phen, float alt_phen)
{
	float ans =  effect_of_interaction + interaction_strength*weighted_effect_of_alternative_individual*(exp(-powf(focal_phen-alt_phen,2.0))/(1+exp(-powf(focal_phen-alt_phen,2.0))));
	
	return ans;
}

inline __host__ __device__ float calculate_effect_of_interaction_noPhen(float effect_of_interaction, float weighted_effect_of_alternative_individual, float interaction_strength)
{
	float ans =  effect_of_interaction + interaction_strength*weighted_effect_of_alternative_individual;
	return ans;
}

inline __host__ __device__ float update_fecundity_phenotype(float effect_of_interaction, float fecundity_phenotype, float trait_value, float fecundity_trait_trade_off)
{
	float ans = fecundity_phenotype * (effect_of_interaction);
	// Implement the pleiotropy process here
	ans = ans*exp(-fecundity_trait_trade_off*powf(trait_value,2.0));
	return ans; 
}

inline __host__ __device__ float update_survivorship_phenotype(float effect_of_interaction, float survivorship_value)
{
	float ans = survivorship_value * (effect_of_interaction);
	return ans;
}


#endif
