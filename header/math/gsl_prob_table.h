#ifndef GSL_PROB_TABLE_H
#define GSL_PROB_TABLE_H

#include <thrust/host_vector.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class GSL_ProbTable
	{
/* 
* A look-up table to be used for simulating from arbitrary discrete distributions.
*/
	public:
		GSL_ProbTable(thrust::host_vector<float>::iterator prob_begin, thrust::host_vector<float>::iterator prob_end);
		~GSL_ProbTable();

		void draw(const gsl_rng * gen, int &result);
	protected:
		void setup(thrust::host_vector<float>::iterator prob_begin, thrust::host_vector<float>::iterator prob_end);

		thrust::host_vector<float> cumulative_prob;

		gsl_ran_discrete_t *alias_table;

		double *probabilities_double;
		int n;			
	};

#endif
