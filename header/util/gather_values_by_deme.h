#ifndef VALUES_BY_DEME
#define VALUES_BY_DEME

#include <util/thrust_functors.h>

#include <iostream>
#include <stdio.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>

/* 

the function gather_values_by_deme works like this:

Suppose you have a vector A of indices

A = [A1, A2, A3, ..., An]

and a corresponding vector of the demes to which the individual Ai belongs

D = [D1, D2, D3 ,..., Dn]

If each deme has some deme-specific value stored in a vector X

B = [X1, X2, ..., Xm]

where m is the number of demes. Then gather_values_by_deme creates a resulting vector C such that:

C = [X[D1], X[D2], ..., X[Dm]]


*/
void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<int> &deme_specific_value, 
			   thrust::device_vector<int> &ans);

void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<float> &deme_specific_value, 
			   thrust::device_vector<float> &ans);

void gather_values_by_deme(thrust::device_vector<int> &indices,
			   thrust::device_vector<int> &demes,
                           thrust::device_vector<float>::iterator deme_specific_values_begin,
                           thrust::device_vector<float> &ans);

#endif
