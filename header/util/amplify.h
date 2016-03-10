#ifndef AMPLIFY_H
#define AMPLIFY_H


#include <iostream>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>


/* If you have a vector of values
A = [A0, A1, A2, A3, ..., An]

and each element gets repeated

B = [B0, B1, B2, B3, ..., Bn] 

times, produce a third vector C that repeats A0 B0 times, A1 B1 times, etc...

e.g., 

A = [0, 1, 2, 3, 4]

B = [5, 2, 5, 0, 1]

then 

C = [0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 4]



*/


/********************************

This function amplify takes two vectors, amplify_counts and values_to_amplify, and stores the results in a vector amplified_values where values_to_amplify[i] is repeated amplify_counts[i] times.

E.g., 

thrust::device_vector<int> amplify_counts(5);
thrust::device_vector<int> values_to_amplify(5);
thrust::device_vector<int> amplified_values(1);

for (int i=0; i < 5; i++)
	{
	values_to_amplify[i] = i;
	amplify_counts[i] = 4 + i;
	}

amplify(amplify_counts, values_to_amplify, amplified_values);

# amplified_values = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
	

********************************/


void amplify(thrust::device_vector<int> &values_to_amplify,
	     thrust::device_vector<int> &amplify_counts,
	     thrust::device_vector<int> &amplified_values);


void amplify_float(thrust::device_vector<float> &values_to_amplify,
	     	   thrust::device_vector<int> &amplify_counts,
	     	   thrust::device_vector<float> &amplified_values);

void amplify_sequence(thrust::device_vector<int> &amplify_counts,
	     int number_of_elements_in_sequence,
	     thrust::device_vector<int> &amplified_values);


#endif
