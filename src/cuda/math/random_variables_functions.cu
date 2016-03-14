#include <math/random_variables_functions.h>
#include <math/thrust_probabilities.h>

#include <curand.h>
#include <fstream>
#include <iostream>

void prime_random_number_generator(curandGenerator_t gen, int seed)
	{
/*
*
* Not clear why having this be in random_variables_functions does not work when prime_random_number_generator gets called from inds_stochastic's initialization. it is likely an issue of code not linking correctly...
*
*/
	int size = 100;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, seed);

	//curand declarations
	thrust::device_vector<float> rand(size);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(gen, rand_ptr, size); // priming up the random number generator takes some time, get it done early.
	rand.clear();
	}

void draw_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
	{
/*
*
* a wrapper for drawing gaussian random variates - curandGenerateNormal() behaves when interacting with thrust vectors. This routine is slower than curandGenerateNormal(), so waiting for workable thrust generator to come out.
*
*/
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	thrust::device_vector<float> rand2(samples_needed); 
	float *rand_ptr2 = raw_pointer_cast(&rand2[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	curandGenerateUniform(gen, rand_ptr2, samples_needed);

	normal_rv gen_normal(mean, stddev);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),rand2.begin(),random_variates.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed,rand2.begin() + samples_needed,random_variates.begin()  + samples_needed)), gen_normal);
	}

void draw_discrete_gaussian(int samples_needed, float mean, float stddev, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	thrust::device_vector<float> rand2(samples_needed); 
	float *rand_ptr2 = raw_pointer_cast(&rand2[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	curandGenerateUniform(gen, rand_ptr2, samples_needed);

	discrete_normal_rv gen_discrete_normal(mean, stddev);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),rand2.begin(),random_variates.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed,rand2.begin() + samples_needed,random_variates.begin()  + samples_needed)), gen_discrete_normal);
	}




void draw_poisson(int samples_needed, float lambda, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
/*
*
* NB: This is something that apparently CUDA 5 also implements. Eventually we should switch to using cuda's native poisson and gaussian simulator
*
*/
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	
	curandGenerateUniform(gen, rand_ptr1, samples_needed);

	thrust::transform(rand1.begin(), rand1.end(), random_variates.begin(), poisson_rv(lambda));
	}

void draw_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	thrust::device_vector<float> rand2(samples_needed); 
	float *rand_ptr2 = raw_pointer_cast(&rand2[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	curandGenerateUniform(gen, rand_ptr2, samples_needed);

	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),rand2.begin(),Means.begin(),Standard_deviations.begin(), random_variates.begin() )),
		thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed,rand2.begin() + samples_needed, Means.begin() + samples_needed, Standard_deviations.begin() + samples_needed, random_variates.begin() + samples_needed)),normal_rv_different_parameters());		
	}


void draw_gaussian_different_parameters(int samples_needed, float mean_value, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> Means(samples_needed);
	thrust::fill(Means.begin(), Means.end(), mean_value);

	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	thrust::device_vector<float> rand2(samples_needed); 
	float *rand_ptr2 = raw_pointer_cast(&rand2[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	curandGenerateUniform(gen, rand_ptr2, samples_needed);

	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),rand2.begin(),Means.begin(),Standard_deviations.begin(), random_variates.begin() )),
		thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed,rand2.begin() + samples_needed, Means.begin() + samples_needed, Standard_deviations.begin() + samples_needed, random_variates.begin() + samples_needed)),normal_rv_different_parameters());		
	}

void draw_discrete_gaussian_different_parameters(int samples_needed, thrust::device_vector<float> &Means, thrust::device_vector<float> &Standard_deviations, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	thrust::device_vector<float> rand2(samples_needed); 
	float *rand_ptr2 = raw_pointer_cast(&rand2[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	curandGenerateUniform(gen, rand_ptr2, samples_needed);

	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),rand2.begin(),Means.begin(),Standard_deviations.begin(), random_variates.begin() )),
		thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed,rand2.begin() + samples_needed, Means.begin() + samples_needed, Standard_deviations.begin() + samples_needed, random_variates.begin() + samples_needed)),discrete_normal_rv_different_parameters());		
	}

void  draw_poisson_different_parameters(int samples_needed, thrust::device_vector<float> &lambdas, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);

	curandGenerateUniform(gen, rand_ptr1, samples_needed);
	
	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),lambdas.begin(), random_variates.begin() )),
		thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed, lambdas.begin() + samples_needed, random_variates.begin() + samples_needed)), poisson_rv_different_parameters());		
	}


void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	curandGenerateUniform(gen, rand_ptr1, samples_needed);

	thrust::fill(random_variates.begin(), random_variates.end(), 0);

	thrust::transform(rand1.begin(), rand1.end(), random_variates.begin(), bernoulli_rv(probability));
	}

void draw_bernoulli(int samples_needed, float probability, thrust::device_vector<float> &random_variates, curandGenerator_t gen)
	{
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	curandGenerateUniform(gen, rand_ptr1, samples_needed);

	thrust::fill(random_variates.begin(), random_variates.end(), 0);

	thrust::device_vector<int> temporary_random_variates(samples_needed);

	thrust::transform(rand1.begin(), rand1.end(), temporary_random_variates.begin(), bernoulli_rv(probability));

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(temporary_random_variates.begin(), random_variates.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(temporary_random_variates.begin()+samples_needed, random_variates.begin()+samples_needed)),
		 int_to_float());
	}

void draw_bernoulli_different_parameters(int samples_needed, thrust::device_vector<float> &probabilities, thrust::device_vector<int> &random_variates, curandGenerator_t gen)
	{
/*
*
* Draw a sequence of samples_needed Bernoulli random variables where the probability of success is different for each random variate
*
*/
	thrust::device_vector<float> rand1(samples_needed); 
	float *rand_ptr1 = raw_pointer_cast(&rand1[0]);
	curandGenerateUniform(gen, rand_ptr1, samples_needed);

	thrust::fill(random_variates.begin(), random_variates.end(), 0);

	thrust::for_each( thrust::make_zip_iterator(thrust::make_tuple(rand1.begin(),probabilities.begin(), random_variates.begin() )),
			thrust::make_zip_iterator(thrust::make_tuple(rand1.begin() + samples_needed, probabilities.begin() + samples_needed, random_variates.begin() + samples_needed)),
			bernoulli_rv_different_parameters());		
	}
