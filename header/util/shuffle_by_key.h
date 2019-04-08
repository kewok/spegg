#ifndef SHUFFLE_BY_KEY
#define SHUFFLE_BY_KEY

#include <curand.h>
#include <iostream>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <util/thrust_functors.h>


/* 
Shuffle by key assumes that if you have a vector of values

V = [V1, V2, V3, ..., Vn]

and a each value has a corresponding key

K = [K1, K2, K3, ..., Kn]

where potentially K_i = K_j for some subset of K, then you can create a random shuffle among the elements of V_i such that all V_i which have key K_i are randomly arranged. 

*/

void shuffle_by_key(thrust::device_vector<int> &keys, thrust::device_vector<int> &values, curandGenerator_t gen, bool shuffle_keys=false);

void shuffle_by_key(thrust::device_vector<int> &keys, thrust::device_vector<float> &values, curandGenerator_t gen);

void shuffle_by_key(thrust::device_vector<float> &keys, thrust::device_vector<int> &values, curandGenerator_t gen);

void shuffle_by_key(thrust::device_vector<float> &keys, thrust::device_vector<float> &values, curandGenerator_t gen);

#endif
