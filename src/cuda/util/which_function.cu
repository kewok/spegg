#include <util/which_function.h>

/* TODO: Right now, which_something are actually predicated on the approach of removing values for which something is false, e.g., which_equal_to(stencil, answer, value) returns indices in vector answer by testing which elements of stencil are not equal to "value". This is potentially confusing.

*/


void which_equal_to(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_not_equal<int> (value)), answer.end() );
	}

void which_equal_to(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value,
	   int elements_evaluated)
	{
	answer.resize( elements_evaluated );

	thrust::device_vector<int> indices_all( elements_evaluated );
	
	thrust::sequence(indices_all.begin(), indices_all.end());

	thrust::device_vector<int> stencil_local( elements_evaluated );
	thrust::copy(stencil.begin(), stencil.begin() + elements_evaluated, stencil_local.begin());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil_local.begin() , answer.begin(),  unary_not_equal<int> (value)), answer.end() );
	}

void which_equal_to(thrust::device_vector<float> &stencil, thrust::device_vector<int> &answer, float value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_not_equal<float> (value)), answer.end() );
	}

void which_greater_than(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_less_equal<int> (value)), answer.end() );
	}

void which_greater_than(thrust::device_vector<float> &stencil,
	   thrust::device_vector<int> &answer,
	   float value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_less_equal<float> (value)), answer.end() );
	}

// Overloaded versions:
void which_greater_than(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value,
	   int elements_evaluated)
	{
	answer.resize( elements_evaluated );

	thrust::device_vector<int> indices_all( elements_evaluated );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	thrust::device_vector<int> stencil_local( elements_evaluated );
	thrust::copy(stencil.begin(), stencil.begin() + elements_evaluated, stencil_local.begin());

	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil_local.begin() , answer.begin(),  unary_less_equal<int> (value)), answer.end() );
	}


void which_less_than(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_greater_equal<int> (value)), answer.end() );
	}

void which_less_than(thrust::device_vector<float> &stencil,
	   thrust::device_vector<int> &answer,
	   float value)
	{
	answer.resize( stencil.size() );

	thrust::device_vector<int> indices_all( stencil.size() );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil.begin() , answer.begin(),  unary_greater_equal<float> (value)), answer.end() );
	}

// Overloaded versions:
void which_less_than(thrust::device_vector<int> &stencil,
	   thrust::device_vector<int> &answer,
	   int value,
	   int elements_evaluated)
	{
	answer.resize( elements_evaluated );

	thrust::device_vector<int> indices_all( elements_evaluated );
	
	thrust::sequence(indices_all.begin(), indices_all.end());
	
	thrust::device_vector<int> stencil_local( elements_evaluated );
	thrust::copy(stencil.begin(), stencil.begin() + elements_evaluated, stencil_local.begin());

	answer.erase(thrust::remove_copy_if(indices_all.begin(), indices_all.end(), stencil_local.begin() , answer.begin(),  unary_greater_equal<int> (value)), answer.end() );
	}

