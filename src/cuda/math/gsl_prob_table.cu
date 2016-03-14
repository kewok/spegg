#include <math/gsl_prob_table.h>

GSL_ProbTable::GSL_ProbTable(thrust::host_vector<float>::iterator prob_begin, thrust::host_vector<float>::iterator prob_end)
	{
	setup(prob_begin, prob_end);
	}

void GSL_ProbTable::setup(thrust::host_vector<float>::iterator prob_begin, thrust::host_vector<float>::iterator prob_end)
	{
	n = thrust::distance(prob_begin, prob_end);
	probabilities_double = (double*) malloc(n*sizeof(double));

	for (int i=0; i < n; i++)
		{
		probabilities_double[i] = (double) *(prob_begin + i);
		}
	alias_table = gsl_ran_discrete_preproc(n, probabilities_double);
	}

void GSL_ProbTable::draw(const gsl_rng * gen, int &result)
	{
	result = gsl_ran_discrete(gen, alias_table);
	}

GSL_ProbTable::~GSL_ProbTable()
	{
	free(probabilities_double);
	gsl_ran_discrete_free(alias_table);
	}
