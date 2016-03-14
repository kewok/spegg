#ifndef UPDATE_FISH_H
#define UPDATE_FISH_H

#include <species/update/updatebehavior.h>
#include "Fish_Habitat.h"

#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

class update_Fish : public UpdateBehavior
	{
	inds_stochastic *species;
	Fish_Habitat *habitat;

	// Constructor
	public:
	update_Fish(inds_stochastic *species, environment *habitat) 
	 	{
		this->species = species;

		this->MORTALITY_PHENOTYPE = species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE"];
		this->FECUNDITY_PHENOTYPE = species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE"];
		this->IRREVERSIBLE_MASS_PHENOTYPE = species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"];
		this->REVERSIBLE_MASS_PHENOTYPE = species->demeParameters->species_specific_values["REVERSIBLE_MASS_PHENOTYPE"];
		this->RESOURCE_LIMITATION_PHENOTYPE = species->demeParameters->species_specific_values["RESOURCE_LIMITATION_PHENOTYPE"];
		this->EGGSIZE_PHENOTYPE =  species->demeParameters->species_specific_values["EGGSIZE_PHENOTYPE"];

		// Copy the constants 
		this->size = species->size;
		this->Number_of_Demes = species->Num_Demes;

		// Nifty hack: Downcast your generic habitat object to a pointer to fish class to access its members
		this->habitat = static_cast<Fish_Habitat*> (habitat);
		this->intra_annual_time_steps = this->habitat->intra_annual_time_steps;
		// Prepare pointers to thrust vectors:
		prepare_growth_constants_pointers();
		prepare_survivorship_constants_pointers();
		}

	void update();

	protected:
		// variable names
		int MORTALITY_PHENOTYPE;
		int FECUNDITY_PHENOTYPE;
		int IRREVERSIBLE_MASS_PHENOTYPE;
		int REVERSIBLE_MASS_PHENOTYPE;
		int RESOURCE_LIMITATION_PHENOTYPE;
		int EGGSIZE_PHENOTYPE;

		int size;
		int Number_of_Demes;
		int intra_annual_time_steps;
		void consume();
		void survive();
		void update_vital_rates();
		void update_fecundity();

		void prepare_growth_constants_pointers();
		void prepare_survivorship_constants_pointers();

		// Pointers to thrust vectors used in survivorship calculation:
		float *effect_of_starvation_ptr;
		float *size_dependent_mortality_constant_ptr;
		float *size_dependent_mortality_coefficient_ptr;

		// Pointers to thrust vectors used in somatic growth calculation:
		float *consumption_allometric_scalar;
		float *gamma;
		float *ontogenetic_niche_shift_constant;
		float *ontogenetic_niche_shift_coefficient;
		float *length_weight_conversion_coefficient;
		float *length_weight_conversion_exponent;

		float *resource_1_maximum;
		float *resource_2_maximum;
		float *handling_time_resource1;
		float *handling_time_resource2;

		float *functional_response_scalar_resource1;
		float *functional_response_scalar_resource2;
		float *alpha_g;

		float *mature_maximum_condition;
		float *juvenile_maximum_condition;

		float *M_sizes_at_maturity;
		float *F_sizes_at_maturity;
	};


/* Functors related to updating */

struct grand_consumption_functor
{
	float *consumption_allometric_scalar;
	float *gamma;
	float *ontogenetic_niche_shift_constant;
	float *ontogenetic_niche_shift_coefficient;
	float *length_weight_conversion_coefficient;
	float *length_weight_conversion_exponent;
	
	float *resource_1;
	float *resource_1_maximum;
	float *resource_2;
	float *resource_2_maximum;
	float *handling_time_resource1;
	float *handling_time_resource2;

	float *functional_response_scalar_resource1;
	float *functional_response_scalar_resource2;

	float *alpha_g;

	float *mature_maximum_condition;
	float *juvenile_maximum_condition;

	float *M_sizes_at_maturity;
	float *F_sizes_at_maturity;
	
	grand_consumption_functor(float *eCoef, float *gammaG, float *ONS_constant, float *ONS_coef, float *lwConv, float *lwExp, float *r1, float *r1_K, float *r2, float *r2_K, float *handling_1, float *handling_2, float *fnl_scalar_1, float *fnl_scalar_2, float *alpha, float *mature_kappa, float *juv_kappa, float *M_size_mat, float *F_size_mat) : consumption_allometric_scalar(eCoef), gamma(gammaG), ontogenetic_niche_shift_constant(ONS_constant), ontogenetic_niche_shift_coefficient(ONS_coef), length_weight_conversion_coefficient(lwConv), length_weight_conversion_exponent(lwExp), resource_1(r1), resource_2(r2), resource_1_maximum(r1_K), resource_2_maximum(r2_K), handling_time_resource1(handling_1), handling_time_resource2(handling_2), functional_response_scalar_resource1(fnl_scalar_1), functional_response_scalar_resource2 (fnl_scalar_2), mature_maximum_condition(mature_kappa), juvenile_maximum_condition(juv_kappa), M_sizes_at_maturity(M_size_mat), F_sizes_at_maturity(F_size_mat), alpha_g(alpha)
	{};

	/* 
		Elements in the tuple.
		----------------------
		0: irreversible_mass
		1: reversible_mass
		2: pop
		3: proportion of adult resource consumed (currently dummy variable)
		4: sex
		5: status
		6: resource 1 consumed
		7: resource 2 consumed
		8: random fluctuations 1
	*/ 
	template <typename tuple>
	__device__
	void operator()(tuple t) {
		int individuals_pop = thrust::get<2>(t);
		
		if (thrust::get<5>(t) > 0)
			{
			// If resources were at maximum abundance that month, and the individual specialized exclusively on this prey type, then the intake rate is integrate_0_30 b W^gamma which comes out to, in Hiyama and Kitahara's notation, the integrated weight increase with no metabolic costs :

			float available_resources_1 = resource_1[individuals_pop]/resource_1_maximum[individuals_pop];
			float satiation_resource_1 = (functional_response_scalar_resource1[individuals_pop]*available_resources_1)/(1 + handling_time_resource1[individuals_pop] * available_resources_1); 

			float available_resources_2 = resource_2[individuals_pop]/resource_2_maximum[individuals_pop];
			float satiation_resource_2 = (functional_response_scalar_resource2[individuals_pop]*available_resources_2)/(1 + handling_time_resource2[individuals_pop] * available_resources_2); 

			float xt0 = thrust::get<0>(t);
			float yt0 = thrust::get<1>(t);

			float xt1 = 0;
			float yt1 = 0;

			float kappa = 0;
			float gamma_val = gamma[individuals_pop];

			float maturity_threshold = 0;

			if (thrust::get<4>(t) == 0)
				maturity_threshold = F_sizes_at_maturity[individuals_pop];
			else
				maturity_threshold = M_sizes_at_maturity[individuals_pop];

			// Set the initial kappa value. 
			// NB: Be very careful with this: "if (x) do_this if(!x) do_that" seems to give different instructions from "if (x) do_this else do_that"
			if (xt0 >= maturity_threshold)
				kappa = 1/((1+mature_maximum_condition[individuals_pop])*mature_maximum_condition[individuals_pop]);
			else
				kappa = 1/((1+juvenile_maximum_condition[individuals_pop])*juvenile_maximum_condition[individuals_pop]);

			//printf("kappa: %f\n", kappa);

			float maximum_consumption = 0;
			float percent_adult_resource_eaten = 0;
			float percent_juvenile_resource = 0;
			float Eg = 0;

			float ons_constant = ontogenetic_niche_shift_constant[individuals_pop];
			float ons_coefficient = ontogenetic_niche_shift_coefficient[individuals_pop];
			float conversion_efficiency = consumption_allometric_scalar[individuals_pop];
			float lw_coefficient = length_weight_conversion_coefficient[individuals_pop];
			float lw_exponent = length_weight_conversion_exponent[individuals_pop];
			
			float alpha_val = alpha_g[individuals_pop];

			// Do the dynamic energy budgeting

			float prey_eaten_1 = 0;
			float prey_eaten_2 = 0;


			maximum_consumption = conversion_efficiency*pow(xt0+yt0, gamma_val);

			//printf("Maximum consumption: %f\n", maximum_consumption);
	// Energy growth is therefore

			float bodyLength = lw_coefficient* pow(xt0,lw_exponent);
				
			percent_adult_resource_eaten = 1/(1+exp(-(ons_constant + ons_coefficient*(bodyLength))));

			percent_juvenile_resource = 1-percent_adult_resource_eaten;
			
			prey_eaten_1 += satiation_resource_1*percent_juvenile_resource*maximum_consumption;
			prey_eaten_2 += satiation_resource_2*percent_adult_resource_eaten*maximum_consumption;

			//printf("Prey eaten juv. v. adult: %f %f ONS juv v. ad.: %f %f Size: %f\n", prey_eaten_1, prey_eaten_2, percent_juvenile_resource, percent_adult_resource_eaten, xt0);

				// Subtract metabolic costs
			Eg = (satiation_resource_1*percent_juvenile_resource + satiation_resource_2*percent_adult_resource_eaten)*maximum_consumption - alpha_val*(xt0+yt0);
	
			if (Eg >= 0)
				{
				xt1 = (yt0/xt0)*kappa*Eg + xt0;
				yt1 = (1-(yt0/xt0)*kappa)*Eg + yt0;
				}
			else
				{
				yt1 = yt0 + Eg;
				if (yt1 < 0)
					{
					yt1 = 0;
					}
				xt1 = xt0;
				}

			//printf("Juvenile Prey available: %f Energy gained: %f Original size: %f %f Final size: %f %f\n", available_resources_1, Eg, xt0, yt0, xt1, yt1);

			thrust::get<0>(t) = xt1;
			thrust::get<1>(t) = yt1;
					
			thrust::get<6>(t) = prey_eaten_1;
			thrust::get<7>(t) = prey_eaten_2;
			}
		}
};


struct calculate_mortality
{
	float *effect_of_starvation;
	float *size_dependent_mortality_constant;
	float *size_dependent_mortality_coefficient;

	calculate_mortality(float *bSTV, float *bMort0, float *bMort1) : effect_of_starvation(bSTV), size_dependent_mortality_constant(bMort0), size_dependent_mortality_coefficient(bMort1)
	{};

	/* 
		Elements in the tuple.


		----------------------
		0: irreversible mass
		1: reversible mass
		2: probability of survivorship
		3: individual's subpopulation
		4: individual's age
		5: individual's dead/alive status
		6: uniform random variable

	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int individuals_pop = thrust::get<3>(t);

		float condition = thrust::get<1>(t)/thrust::get<0>(t);

		if (thrust::get<5>(t)==1) /* If the individual is alive */
			{
			if (thrust::get<4>(t) >= 0) // If individual is not an egg
				{
				//thrust::get<3>(t) =(1-exp(-effect_of_starvation[individuals_pop]*condition))* (1-exp(size_dependent_mortality_coefficient[individuals_pop]*sqrt(thrust::get<1>(t)) + size_dependent_mortality_constant[individuals_pop]));
				// Alternative formulation based on logistic regression is:
				thrust::get<2>(t) =(1-exp(-effect_of_starvation[individuals_pop]*condition))* (1/(1+exp(size_dependent_mortality_coefficient[individuals_pop]*(thrust::get<0>(t) + size_dependent_mortality_constant[individuals_pop]))));
				}

			if (thrust::get<2>(t) > 1)
				{
				thrust::get<2>(t) = 1;
				}

			thrust::get<2>(t) = thrust::get<2>(t) * 0.9995; // Impose a maximum background survivorship -> leads to mean life span of ~25 yrs = 1/(25*90)
			
			if (thrust::get<2>(t) < 0)	
				{
				thrust::get<2>(t)=0;
				}
			/*
			float ans;
			ans = thrust::get<2>(t);
			float size = thrust::get<0>(t);

			if (size < 10)
				printf("Condition: %f Body size: %f Reversible_mass: %f Daily Surivovrship: %f\n", condition, size, thrust::get<1>(t), ans );
			*/

			// Actually kill them -> possibly have a different functor?
			if (thrust::get<6>(t) < thrust::get<2>(t))
				thrust::get<5>(t) = 1;
			else
				thrust::get<5>(t) = 0;
		
			}
		}	
};


// Update fecundity according to (reversible-irreversible*max_cond)/(max_cond*eggsize + egg_size)

struct update_female_fecundity_functor
{
	float *maximum_condition;
	
	update_female_fecundity_functor(float *maximum_conditions_ptr) : maximum_condition(maximum_conditions_ptr)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: sex of the individual
		1: the individual's subpopulation
		2: the individual's irreversible mass
		3: the individual's reversible mass
		4: the genetically determined average offspring size
		5: the individual's final fecundity score
		6: the individual's index
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float eggsize = thrust::get<4>(t);
		thrust::get<5>(t) = 0; // make sure there is no carry-over from the previous time step.
		thrust::get<5>(t) = (1-thrust::get<0>(t)) * (thrust::get<3>(t) - thrust::get<2>(t)*maximum_condition[thrust::get<1>(t)])/(maximum_condition[thrust::get<1>(t)]*eggsize + eggsize);
		// because we are dealing with the number of kids, convert to integers
		int tempVal = (int) thrust::get<5>(t);
		thrust::get<5>(t) = (float) tempVal;

		float surplus_mass_used_in_reproduction = (1-thrust::get<0>(t)) *(thrust::get<3>(t) - thrust::get<2>(t)*maximum_condition[thrust::get<1>(t)]);


		float finalfecund = thrust::get<5>(t);
		float irrevmass = thrust::get<2>(t);
		float revmass = thrust::get<3>(t);
		float offsize = thrust::get<4>(t);
		/*
		if (finalfecund > 0)
			printf("%f %f %f %f %f\n",finalfecund, irrevmass, revmass, surplus_mass_used_in_reproduction, offsize);
		*/
		/* make sure that you deduct the amount used from reproductive females */
		if (thrust::get<5>(t) > 0)
			thrust::get<3>(t) -= surplus_mass_used_in_reproduction;
		else
			thrust::get<5>(t) = 0 ;
		}
};

#endif
