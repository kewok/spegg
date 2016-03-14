#ifndef RANDOM_VARIABLES_H
#define RANDOM_VARIABLES_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <util/thrust_functors.h>

void prime_random_number_generator(curandGenerator_t gen, int seed);

void draw_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<float> &random_variates, curandGenerator_t gen);

void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

void draw_poisson(int samples_needed, float lambda, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

void draw_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen);

void draw_gaussian_different_parameters(int samples_needed, float mean_value, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen);

void draw_discrete_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

void  draw_poisson_different_parameters(int samples_needed, thrust::device_vector<float> &lambdas, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<float> &random_variates, curandGenerator_t gen);

void draw_bernoulli_different_parameters(int samples_needed, thrust::device_vector<float> &probabilities, thrust::device_vector<int> &random_variates, curandGenerator_t gen);

#endif

