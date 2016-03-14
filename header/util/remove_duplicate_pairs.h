#ifndef REMOVE_DUPLICATE_PAIRS_H
#define REMOVE_DUPLICATE_PAIRS_H

#include <thrust/host_vector.h>

/*************************
//
//
// util function for identifying and removing duplicates in pairs made up of two vectors; e.g.
// vectorA = [0,0,0, 1,1,1, 2,2,2, 0,0]
// vectorB = [3,4,5, 6,6,7, 8,9,9, 6,3]
// 
// Remove pairs so that:
//
// vectorA = [0,0,0,0, 1,1, 2,2]
// vectorB = [3,4,5,6, 6,7, 8,9]
//
// Note that vectorA and vectorB must be of the same size.
//
*************************/

void remove_duplicate_pairs(thrust::host_vector<int> &vectorA, 
			    thrust::host_vector<int> &vectorB );
#endif
