#include <math/random_variables_functions.h>
#include <math/thrust_probabilities.h>

#include <fstream>
#include <iostream>

void prime_random_number_generator(gsl_rng* gen, int seed)
	{
/*
*Initialize the random number generator
*/	
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	gen = gsl_rng_alloc (T);
	gsl_rng_set(gen, seed);
	}

void draw_gaussian(int samples_needed, float mean, float stddev, thrust::host_vector<float> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = mean + gsl_ran_gaussian(gen, stddev);
		}
	}

void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::host_vector<int> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (int) (mean + gsl_ran_gaussian(gen, stddev));
		}
	}


void draw_poisson(int samples_needed, float lambda, thrust::host_vector<int> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (int) (gsl_ran_poisson(gen, lambda));
		}
	}

void draw_gaussian_different_parameters(int samples_needed, thrust::host_vector<float> &Means, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<float> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (Means[i] + gsl_ran_gaussian(gen, Standard_deviations[i]));
		}	
	}


void draw_gaussian_different_parameters(int samples_needed, float mean_value, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<float> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (mean_value + gsl_ran_gaussian(gen, Standard_deviations[i]));
		}
	}

void draw_discrete_gaussian_different_parameters(int samples_needed, thrust::host_vector<float> &Means, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<int> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (int) (Means[i] + gsl_ran_gaussian(gen, Standard_deviations[i]));
		}	
	}

void  draw_poisson_different_parameters(int samples_needed, thrust::host_vector<float> &lambdas, thrust::host_vector<int> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (int) (gsl_ran_poisson(gen, lambdas[i]));
		}		
	}


void draw_bernoulli(int samples_needed, float probability, thrust::host_vector<int> &random_variates, gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (int) (gsl_ran_bernoulli(gen, probability));
		}
	}

void draw_bernoulli(int samples_needed, float probability, thrust::host_vector<float> &random_variates,  gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (float) (gsl_ran_bernoulli(gen, probability));
		}
	}

void draw_bernoulli_different_parameters(int samples_needed, thrust::host_vector<float> &probabilities, thrust::host_vector<float> &random_variates,  gsl_rng *  gen)
	{
	for (int i=0; i < random_variates.size(); i++)
		{
		random_variates[i] = (float) (gsl_ran_bernoulli(gen, probabilities[i]));
		}
	}
