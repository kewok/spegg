#ifndef RANDOM_VARIABLES_H
#define RANDOM_VARIABLES_H

#include <thrust/host_vector.h>
#include <util/thrust_functors.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void draw_gaussian(int samples_needed, float mean, float stddev, thrust::host_vector<float> &random_variates, gsl_rng *  gen);

void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::host_vector<int> &random_variates, gsl_rng *  gen);

void draw_poisson(int samples_needed, float lambda, thrust::host_vector<int> &random_variates, gsl_rng *  gen);

void draw_gaussian_different_parameters(int samples_needed, thrust::host_vector<float> &Means, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<float> &random_variates, gsl_rng *  gen);

void draw_gaussian_different_parameters(int samples_needed, float mean_value, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<float> &random_variates, gsl_rng *  gen);

void draw_discrete_gaussian_different_parameters(int samples_needed, thrust::host_vector<float> &Means, thrust::host_vector<float> &Standard_deviations, thrust::host_vector<int> &random_variates, gsl_rng *  gen);

void  draw_poisson_different_parameters(int samples_needed, thrust::host_vector<float> &lambdas, thrust::host_vector<int> &random_variates, gsl_rng *  gen);

void draw_bernoulli(int samples_needed, float probability, thrust::host_vector<int> &random_variates, gsl_rng *  gen);

void draw_bernoulli(int samples_needed, float probability, thrust::host_vector<float> &random_variates, gsl_rng *  gen);

void draw_bernoulli_different_parameters(int samples_needed, thrust::host_vector<float> &probabilities, thrust::host_vector<int> &random_variates, gsl_rng *  gen);
#endif

