#include "prey_class.h"
#include <math/random_variables_functions.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <thrust/functional.h>

typedef thrust::host_vector<float>::iterator floatIter_t;
typedef thrust::host_vector<int>::iterator intIter_t;

prey::prey(int num_demes, int prey_index, int seed_val, int num_time_steps) : Num_Demes(num_demes), Prey_Index(prey_index), seed (seed_val),  Number_of_Timesteps (num_time_steps)
	{
	// Data vectors of length Num_Demes
	prey_abundance.resize(Num_Demes);
	prey_maximum_abundance.resize(Num_Demes);
	prey_assimilation_efficiency.resize(Num_Demes);

	prey_carrying_capacity.resize(Num_Demes);
	prey_unconstrained_growth_rates.resize(Num_Demes);

	prey_growth_rate_noise_stddev.resize(Num_Demes);

	//Initialize random number generator.
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	gen = gsl_rng_alloc (T);
	gsl_rng_set(gen, seed);

	prey_abundance.resize(Num_Demes);

	//Read in deme-specific data	
	preyParameters = new PreyDemeSpecificData("environment_config.txt", prey_index);
	}


void prey::specify_prey_properties_by_deme()
	{
	thrust::copy(preyParameters->get_vector_ptr("prey_carrying_capacity"), preyParameters->get_vector_ptr("prey_carrying_capacity") + Num_Demes, prey_carrying_capacity.begin());

	thrust::copy(preyParameters->get_vector_ptr("prey_unconstrained_growth_rates"), preyParameters->get_vector_ptr("prey_unconstrained_growth_rates") + Num_Demes, prey_unconstrained_growth_rates.begin());

	thrust::copy(preyParameters->get_vector_ptr("prey_assimilation_efficiency"), preyParameters->get_vector_ptr("prey_assimilation_efficiency") + Num_Demes, prey_assimilation_efficiency.begin());

	thrust::copy(preyParameters->get_vector_ptr("prey_growth_rate_noise_stddev"), preyParameters->get_vector_ptr("prey_growth_rate_noise_stddev") + Num_Demes, prey_growth_rate_noise_stddev.begin());

	thrust::copy(preyParameters->get_vector_ptr("prey_maximum_abundance"), preyParameters->get_vector_ptr("prey_maximum_abundance") + Num_Demes, prey_maximum_abundance.begin());

	// initialize the prey abundances.
	thrust::copy(prey_maximum_abundance.begin(), prey_maximum_abundance.begin() + Num_Demes, prey_abundance.begin());
	}

void prey::update_prey_abundance(thrust::host_vector<float> &effect_of_individuals_on_prey)
	{
	//typedef for clarity.
	typedef thrust::tuple<floatIter_t, floatIter_t, floatIter_t, floatIter_t> tuple_t;
	typedef thrust::zip_iterator<tuple_t> zipIter_t;

	tuple_t start = thrust::make_tuple(prey_unconstrained_growth_rates.begin(), prey_carrying_capacity.begin(), effect_of_individuals_on_prey.begin(), prey_abundance.begin());
	tuple_t end = thrust::make_tuple(prey_unconstrained_growth_rates.begin()  + Num_Demes, prey_carrying_capacity.begin()  + Num_Demes, effect_of_individuals_on_prey.begin() + Num_Demes, prey_abundance.begin() + Num_Demes);

 	//Create zip iterators.
	zipIter_t zstart = thrust::make_zip_iterator(start);
	zipIter_t zend = thrust::make_zip_iterator(end);

	thrust::for_each(zstart, zend,  update_prey());

	// Add random noise to the prey abundance
	thrust::host_vector<float> stochastic_component(Num_Demes);
	thrust::host_vector <float> zeros(Num_Demes);
	thrust::fill(zeros.begin(), zeros.begin() + Num_Demes, 0);

	for (int i=0; i < Num_Demes; i++)
		{
		stochastic_component[i] = gsl_ran_gaussian(gen, prey_growth_rate_noise_stddev[i]);
		}

	thrust::transform(prey_abundance.begin(), prey_abundance.begin() + Num_Demes, stochastic_component.begin(), prey_abundance.begin(), thrust::plus<float>());	

	thrust::replace_if(prey_abundance.begin(), prey_abundance.begin() + Num_Demes, is_less_than_zero_f(), 0);
	}

