#ifndef REDUCE_BY_KEY_WITH_ZEROS_H
#define REDUCE_BY_KEY_WITH_ZEROS_H

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/unique.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/extrema.h>

void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<int> &values, thrust::device_vector<int> &values_output, int number_of_values, int number_of_possible_keys);

void reduce_by_key_with_zeros(thrust::device_vector<int> &key_vector, thrust::device_vector<float> &values, thrust::device_vector<float> &values_output, int number_of_values, int number_of_possible_keys);

#endif
