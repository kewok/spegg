#ifndef INTERACTION_KERNEL_FUNCTORS_H
#define INTERACTION_KERNEL_FUNCTORS_H

#include <curand.h>
#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/random.h>

#define PI 3.14159265358979f
#define POISSON_MAX_RV 0.99999991f
#define MAXIMUM_ITERATIONS_REJECTION_SAMPLING 1000000
#define MAX_SAMPLED 20

inline __host__ __device__ int draw_poisson_rv(float lambda, float use_rv)
{
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


inline __host__ __device__ int discrete_normal_rv(float lambda, float use_rv, float use_rv2)
{
	float r = sqrtf(-2.0f * logf(use_rv));
	float phi = 2 * PI * use_rv2;
	int ans = (int) sqrtf(lambda) * r * cosf(phi) + lambda;
	return ans;	
}

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

/* Functors implementing the interactions */

struct interaction_kernel
{
float *focal_phen1;
float *alt_phen1;
int *alt_status;
int *cumul_deme_size;

int species1_IDX;
int species2_IDX;

interaction_kernel(int focal_species_ID, int alternative_species_ID, float *focal_phen1_ptr, float *alt_phen1_ptr, int *alt_status_ptr, int *cumul_deme_size_ptr) : species1_IDX(focal_species_ID), species2_IDX(alternative_species_ID), focal_phen1(focal_phen1_ptr), alt_phen1(alt_phen1_ptr), alt_status(alt_status_ptr), cumul_deme_size(cumul_deme_size_ptr)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: individual index
		1: individuals deme
		2: average number of other individuals sampled
		3: focal species individuals' status
		4: focal species individuals' fecundity
		5: focal species individuals' survivorship
		6: deme-wide interaction strength for the effect on fecundity
		7: deme-wide interaction strength for the effect on survivorship
		8: trade-off for trait-value versus fecundity
		9: random number seed
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<3>(t) > 0) // If the individual is alive
			{				
			unsigned int seedval = thrust::get<9>(t);
			thrust::default_random_engine rng(seedval);
			thrust::uniform_real_distribution<float> u01(0.0,1.0);

			float use_rv = u01(rng);

			int number_others_sampled = 0;
	
			float lambda = thrust::get<2>(t);
			if (lambda > 20.0)
				{
				float use_rv2 = u01(rng);
				number_others_sampled = discrete_normal_rv(lambda, use_rv, use_rv2);
				}
			else
				{
				number_others_sampled = draw_poisson_rv(lambda, use_rv);
				}
			if (number_others_sampled < 0)
				{
				number_others_sampled = 0;
				}

			int focal_individual = thrust::get<0>(t);
			int ind_deme = thrust::get<1>(t);
			//printf("The number of other individuals sampled by %d is: %d\t in deme %d\n", focal_individual, number_others_sampled, ind_deme);

			int min_ind_index, max_ind_index;
			if (thrust::get<1>(t)==0)
				{
				min_ind_index = 0;
				}
			else
				{
				min_ind_index = cumul_deme_size[ind_deme-1];
				}
			max_ind_index = cumul_deme_size[ind_deme];

			// set the next two parameters to 0 if additive
			float effect_of_interaction_on_fecundity = 1; 
			float effect_of_interaction_on_survivorship = 1;

			float interaction_strength_for_fecundity = thrust::get<6>(t);
			float interaction_strength_for_survivorship = thrust::get<7>(t);
			
			// Actual calculations
			if (min_ind_index != max_ind_index) // If there is at least one other individual that is alive 
				{
				for (int i=0; i < number_others_sampled; i++)
					{
					use_rv = u01(rng);
					int alternative_individual = draw_ind(min_ind_index, max_ind_index, use_rv);
					//printf("The individual encountered by %d is: %d with current effect size %f\n", focal_individual, alternative_individual, effect_of_interaction_on_fecundity);
				
					float weighted_effect_of_alternative_individual = 1;


					
					// Discount effect if the individual from the alternative species is dead. Note that alternative suitability criteria can be used here; e.g., if the alternative individual is a male, if the alternative individual is yourself, etc...
					if ( (alt_status[alternative_individual] == 0) || ((species1_IDX == species2_IDX) && (alternative_individual == focal_individual) ) )
						{
						weighted_effect_of_alternative_individual = 0;	
						}

					// For heterospecifics, the interaction strengths are mediated by trait values focal_phen1 and alt_phen1
					if ((species1_IDX != species2_IDX))
						{
						effect_of_interaction_on_fecundity = calculate_effect_of_interaction(effect_of_interaction_on_fecundity, weighted_effect_of_alternative_individual, interaction_strength_for_fecundity, focal_phen1[focal_individual], alt_phen1[alternative_individual]);
						effect_of_interaction_on_survivorship = calculate_effect_of_interaction(effect_of_interaction_on_survivorship, weighted_effect_of_alternative_individual, interaction_strength_for_survivorship, focal_phen1[focal_individual], alt_phen1[alternative_individual]);
						//printf("Phenotypes involved between %d and %d are: %f %f with effect %f and weighted effect %f\n", focal_individual, alternative_individual, focal_phen1[focal_individual], alt_phen1[alternative_individual], effect_of_interaction_on_fecundity, weighted_effect_of_alternative_individual);
						}

					// For conspecific, the interactions strengths are not mediated by the traits
					else
						{
						effect_of_interaction_on_fecundity = calculate_effect_of_interaction_noPhen(effect_of_interaction_on_fecundity, weighted_effect_of_alternative_individual, interaction_strength_for_fecundity);
						effect_of_interaction_on_survivorship = calculate_effect_of_interaction_noPhen(effect_of_interaction_on_survivorship, weighted_effect_of_alternative_individual, interaction_strength_for_survivorship);
						}
					}
				}
						
			// clean up the random number generator
			rng.discard(seedval);

			float old_fecundity = thrust::get<4>(t);
			float fecundity_trait_trade_off = thrust::get<8>(t);

			//printf("The outcome for individual %d is: %f %f %f %f\n", focal_individual, effect_of_interaction_on_fecundity, old_fecundity, fecundity_trait_trade_off, interaction_strength_for_fecundity);

			// Verify in R via:
			// fu <- function(interaction_strength, weighted_effect_of_alternative_individual, focal_phen, alt_phen) interaction_strength + interaction_strength*weighted_effect_of_alternative_individual*(exp(-powf(focal_phen-alt_phen,2.0))/(1+exp(-powf(focal_phen-alt_phen,2.0)))); note there might be some minor discrepency due to using floating point precision.

			thrust::get<4>(t) = update_fecundity_phenotype(effect_of_interaction_on_fecundity, old_fecundity, focal_phen1[focal_individual], fecundity_trait_trade_off);

			float new_fecundity = thrust::get<4>(t);
			//printf("The new fecundity of individual %d is: %f\t", focal_individual, new_fecundity);

			float old_survivorship = thrust::get<5>(t);

			thrust::get<5>(t) = update_survivorship_phenotype(effect_of_interaction_on_survivorship, old_survivorship);

			float new_survivorship = thrust::get<5>(t);
			//printf("Their new survivorship of individual %d is: %f\n", focal_individual, new_survivorship);
			}
		}
	};

#endif
