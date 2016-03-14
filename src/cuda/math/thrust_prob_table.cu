#include <math/thrust_prob_table.h>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>

void ThrustProbTable::setup(thrust::device_vector<float>::iterator prob_begin, thrust::device_vector<float>::iterator prob_end)
	{
	int n = thrust::distance(prob_begin, prob_end);
	cumulative_prob.resize(n);
	
	float total = thrust::reduce(prob_begin, prob_end);
	thrust::device_vector<float> total_vec(n);
	thrust::fill(total_vec.begin(), total_vec.end(), total);
	
	thrust::device_vector<float> relative_prob(n);
	thrust::transform(prob_begin, prob_end, total_vec.begin(), relative_prob.begin(), thrust::divides<float>());
	
	thrust::inclusive_scan(relative_prob.begin(), relative_prob.end(), cumulative_prob.begin());
	}

void ThrustProbTable::draw(thrust::device_vector<float>::iterator uniform_begin, thrust::device_vector<float>::iterator uniform_end, thrust::device_vector<int>::iterator result)
	{
	thrust::lower_bound(cumulative_prob.begin(), cumulative_prob.end(), uniform_begin, uniform_end, result);
	}

thrust::device_ptr<float> ThrustProbTable::get_cumulative_probabilities_ptr()
	{
	return(&cumulative_prob[0]);
	}
