#include <util/one_dim_two_dim.h>

struct one_dim_two_dim_functor
	{
/* 
*For a one dimensional vector made by converting a two-dimensional * vector into a 1D vector, find the corresponding values from the original 2-D vector. 
If you have a vector of int values
A = [A0, A1, A2, A3, ..., An]

and another vector of int values 

B = [B0, B1, B2, B3, ..., Bm] 

then form a third vector

C = [C0, C1, ...., Cm, C(m+1), ...,C(nxm)]

where C0 corresponds to [A0,B0], C1 corresponds to [A0,B1], C2 to [A0,B2], ...,Cm to [A0,Bm], C(m+1) to [A1,B0],..., C(nxm) to [An,Bm].

Conversely, if you have a vector

C = [C0, C1, ...., Cm, C(m+1), ...,C(nxm)]

then one_dim_two_dim works such that for each C_j, return the pair [Ai,Bk] that correspond to C_j, e.g., one_dim_two_dim(C(m+1), length(A), length(B))=[A1,B0].

Note that all of these are o(1) operations per variable element.
*/

	int length_vector_2;

	one_dim_two_dim_functor(int len_vector_2) : length_vector_2(len_vector_2)
	{};
	/* 
		Elements in the tuple.

		----------------------
		0: input value from the single-dimensioned array
		1: corresponding output value for the first array
		2: corresponding output value for the second array
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float one_dim_val = (float) thrust::get<0>(t);
		float len_vec_2 = (float) length_vector_2;
		thrust::get<1>(t) = floorf(one_dim_val/len_vec_2);
		thrust::get<2>(t) = thrust::get<0>(t) - thrust::get<1>(t)*length_vector_2;
		}
	};

void one_dim_two_dim(thrust::device_vector<int> &vector1_values,
	     thrust::device_vector<int> &vector2_values,
	     thrust::device_vector<int> &new_vector,
	     thrust::device_vector<int> &values_for_vector_1,
	     thrust::device_vector<int> &values_for_vector_2)
	{

	/* If you have a vector of int values
	*A = [A0, A1, A2, A3, ..., An]
	*and another vector of int values 
	*B = [B0, B1, B2, B3, ..., Bm] 
	*then form a third vector
	*C = [C0, C1, ...., Cm, C(m+1), ...,C(nxm)]
	*where C0 corresponds to [A0,B0], C1 corresponds to [A0,B1], C2 to [A0,B2], ...,Cm to [A0,Bm], C(m+1) to [A1,B0],..., C(nxm) to [An,Bm].
	*Conversely, if you have a vector
	*C = [C0, C1, ...., Cm, C(m+1), ...,C(nxm)]
	*then one_dim_two_dim works such that for each C_j, return the pair [Ai,Bk] that correspond to C_j, e.g., one_dim_two_dim(C(m+1), length(A), length(B))=[A1,B0].
	*Note that all of these are o(1) operations per variable element.
	*/
	int length_new_vector = new_vector.size();
	int length_vector_2 = vector2_values.size();

	thrust::device_vector<int> indices_vector_1 ( length_new_vector );
	thrust::device_vector<int> indices_vector_2 ( length_new_vector );

	one_dim_two_dim_functor get_two_dim(length_vector_2);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(new_vector.begin(), indices_vector_1.begin(), indices_vector_2.begin())), 	
			 thrust::make_zip_iterator(thrust::make_tuple(new_vector.end(), indices_vector_1.end(), indices_vector_2.end())),
			 get_two_dim);

	thrust::gather(indices_vector_1.begin(), indices_vector_1.end(), vector1_values.begin() , values_for_vector_1.begin());
	thrust::gather(indices_vector_2.begin(), indices_vector_2.end(), vector2_values.begin() , values_for_vector_2.begin());
	}

