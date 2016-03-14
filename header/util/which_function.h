#ifndef WHICH_FUNCTION_H
#define WHICH_FUNCTION_H

#include <util/thrust_functors.h>

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>

// give the ‘TRUE’ indices of a stencil, creating array indices.

void which_equal_to(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value);

void which_equal_to(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value);

void which_greater_than(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value);

void which_greater_than(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value);

void which_less_than(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value);

void which_less_than(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value);


// overloaded versions, which consider only elements up to elements_evaluated 
void which_equal_to(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value, int elements_evaluated);

void which_equal_to(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value, int elements_evaluated);

void which_greater_than(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value, int elements_evaluated);

void which_greater_than(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value, int elements_evaluated);

void which_less_than(thrust::device_vector<int> &stencil, thrust::device_vector<int> &answer, int value, int elements_evaluated);

void which_less_than(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value, int elements_evaluated);

#endif
