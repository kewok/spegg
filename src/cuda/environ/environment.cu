#include <environ/environment.h>

#include <fstream>
#include <iostream>
#include <thrust/count.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/scatter.h>

environment::environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes) : seed(seed_val), ndemes(num_demes), nbiotic_vars(num_biotic_variables), nabiotic_vars(num_abiotic_variables)
	{
	//Allocate biotic data vectors.
	biotic_variables = new thrust::host_vector<float>[nbiotic_vars];
	effect_of_inds_on_biotic_variable = new thrust::host_vector<float>[nbiotic_vars];

	// Resize to accord with number of demes
	for (int i = 0 ; i < nbiotic_vars ; i++) 
		{
		biotic_variables[i].resize(ndemes);
		effect_of_inds_on_biotic_variable[i].resize(ndemes);
		}

		// Initialize the first part of the feedback process 
	for (int i = 0 ; i < nbiotic_vars ; i++) 
		{
		thrust::fill(effect_of_inds_on_biotic_variable[i].begin(), effect_of_inds_on_biotic_variable[i].begin() + ndemes, 0);
		}

	// Allocate abiotic data vectors
	abiotic_variables = new thrust::host_vector<float>[nabiotic_vars];
	
	for (int i=0; i < nabiotic_vars; i++)
		{
		abiotic_variables[i].resize(ndemes);
		}

	//Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	gen = gsl_rng_alloc (T);
	gsl_rng_set(gen, seed);
	}

float* environment::get_abiotic_vector_ptr(const char *abiotic_variable_name)
	{
	return(&abiotic_variables[abiotic_variable_indices[abiotic_variable_name]][0]);
	}


environment::~environment()
	{
	delete[] biotic_variables;
	delete[] effect_of_inds_on_biotic_variable;
	}


