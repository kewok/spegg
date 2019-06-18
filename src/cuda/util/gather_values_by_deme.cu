#include <util/gather_values_by_deme.h>
#include <util/amplify.h>

/* 

the function gather_values_by_deme works like this:

Suppose you have a vector A of indices

A = [A1, A2, A3, ..., An]

and a corresponding vector of any of the m demes (where m is the number of demes) to which the individual Ai belongs

D = [D1, D2, D3 ,..., Dn]

If each deme has some deme-specific value stored in a vector X

B = [X1, X2, ..., Xm].

Then gather_values_by_deme creates a resulting vector C such that:

C = [X[D1], X[D2], ..., X[Dm]]


*/


void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<int> &deme_specific_value, 
			   thrust::device_vector<int> &ans)
	{
	thrust::device_vector<int> demes_of_indices(indices.size());
	thrust::gather(indices.begin(), indices.end(), demes.begin(), demes_of_indices.begin());

	ans.resize(indices.size());
	
	thrust::gather(demes_of_indices.begin(), demes_of_indices.end(), deme_specific_value.begin(), ans.begin());	
	}


void gather_values_by_deme(thrust::device_vector<int> &indices, 
			   thrust::device_vector<int> &demes, 
			   thrust::device_vector<float> &deme_specific_value, 
			   thrust::device_vector<float> &ans)
	{
	thrust::device_vector<int> demes_of_indices(indices.size());
	thrust::gather(indices.begin(), indices.end(), demes.begin(), demes_of_indices.begin());

	ans.resize(indices.size());
	
	thrust::gather(demes_of_indices.begin(), demes_of_indices.end(), deme_specific_value.begin(), ans.begin());	
	}

void gather_values_by_deme(thrust::device_vector<int> &indices,
			   thrust::device_vector<int> &demes,
			   thrust::device_vector<float>::iterator deme_specific_values_begin,
			   thrust::device_vector<float> &ans)
	{
	thrust::device_vector<int> demes_of_indices(indices.size());
	thrust::gather(indices.begin(), indices.end(), demes.begin(), demes_of_indices.begin());

	ans.resize(indices.size());

	thrust::gather(demes_of_indices.begin(), demes_of_indices.end(), deme_specific_values_begin, ans.begin());
	}
