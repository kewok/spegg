#ifndef ONE_DIM_TWO_DIM_H
#define ONE_DIM_TWO_DIM_H

#include <iostream>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>
#include <thrust/functional.h>

void one_dim_two_dim(thrust::host_vector<int> &vector1_values,
	     thrust::host_vector<int> &vector2_values,
	     thrust::host_vector<int> &new_vector,
	     thrust::host_vector<int> &values_for_vector_1,
	     thrust::host_vector<int> &values_for_vector_2);

/* one_dim_two_dim can be tested in R via:
A <- 0:19; B <- 0:29
foo <- function(x)
{
ans1 <- floor(x/length(B))
ans2 <- x-ans1*length(B)
return(c(ans1,ans2))
}
foo(438)
*/

#endif
