#ifndef ASSORTATIVE_MATING_NEONATES_H
#define ASSORTATIVE_MATING_NEONATES_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include "footimer2.h"
#include "inds.h"
#include "reduce_by_key_with_zeroes.h"
#include "deme_specific_data_class.h"
#include "amplify.h"
#include "mating_subpop_thrust_prob_table.h"
#include "random_variables_functions.h"
#include "neonates_class.h"
#include "one_dim_two_dim.h"

#include <thrust/adjacent_difference.h>

class Assortative_mating_neonates :  public EggsNeonates 
{
public:
	Assortative_mating_neonates(thrust::device_vector<int> &pair_populations, DemeSettings *subpopParameters,
			   thrust::device_vector<int> &everybodys_deme,
			   thrust::device_vector<int> &kids_per_mom,
			   thrust::device_vector<int> &current_deme_sizes,
			   thrust::device_vector<int> &maximum_deme_sizes,
			   int N_alive_inds,
			   int num_loci,
			   int nPhen);

	void inherit_genotypes_by_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
				thrust::device_vector<int> &fathers_list,
				thrust::device_vector<int> &mothers_list,
				thrust::device_vector<float> *&fgenotype,
				thrust::device_vector<float> *&mgenotype,
				curandGenerator_t generator);

	void get_mating_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
			      thrust::device_vector<int> &fathers_list,
			      thrust::device_vector<int> &mothers_list,
			      curandGenerator_t generator);


	void record_parents(thrust::device_vector<int> &maternal_id, 
			     thrust::device_vector<int> &paternal_id,
			     thrust::device_vector<int> &ids);

protected:
	thrust::device_vector<int> mothers_chosen;
	thrust::device_vector<int> fathers_chosen;
 

	void get_maternally_derived_genotype_deterministic(thrust::device_vector<int> &mother_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype,
					     curandGenerator_t generator);

	void get_paternally_derived_genotype_deterministic(thrust::device_vector<int> &father_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype,
					     curandGenerator_t generator);

	thrust::device_vector<int> pairs_per_deme;
};
#endif
