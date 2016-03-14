#ifndef ASSORTATIVE_MATING_NEONATES_H
#define ASSORTATIVE_MATING_NEONATES_H

#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <util/footimer2.h>
#include <species/inds.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <species/deme_specific_data_class.h>
#include <util/amplify.h>
#include <math/mating_subpop_thrust_prob_table.h>
#include <math/random_variables_functions.h>
#include <species/add_kids/neonates_class.h>
#include <util/one_dim_two_dim.h>

#include <thrust/adjacent_difference.h>

class Assortative_mating_neonates :  public EggsNeonates 
{
public:
	Assortative_mating_neonates(thrust::host_vector<int> &pair_populations, DemeSettings *subpopParameters,
			   thrust::host_vector<int> &everybodys_deme,
			   thrust::host_vector<int> &kids_per_mom,
			   thrust::host_vector<int> &current_deme_sizes,
			   thrust::host_vector<int> &maximum_deme_sizes,
			   int N_alive_inds,
			   int num_loci,
			   int nPhen);

	void inherit_genotypes_by_pair(thrust::host_vector<float> &probability_pair_becomes_parents,
				thrust::host_vector<int> &fathers_list,
				thrust::host_vector<int> &mothers_list,
				thrust::host_vector<float> *&fgenotype,
				thrust::host_vector<float> *&mgenotype,
				gsl_rng* generator);

	void get_mating_pair(thrust::host_vector<float> &probability_pair_becomes_parents,
			      thrust::host_vector<int> &fathers_list,
			      thrust::host_vector<int> &mothers_list,
			      gsl_rng* generator);


	void record_parents(thrust::host_vector<int> &maternal_id, 
			     thrust::host_vector<int> &paternal_id,
			     thrust::host_vector<int> &ids);

protected:
	thrust::host_vector<int> mothers_chosen;
	thrust::host_vector<int> fathers_chosen;
 

	void get_maternally_derived_genotype_deterministic(thrust::host_vector<int> &mother_index,
					     thrust::host_vector<float> *&mgenotype,
					     thrust::host_vector<float> *&fgenotype,
					     gsl_rng* generator);

	void get_paternally_derived_genotype_deterministic(thrust::host_vector<int> &father_index,
					     thrust::host_vector<float> *&mgenotype,
					     thrust::host_vector<float> *&fgenotype,
					     gsl_rng* generator);

	thrust::host_vector<int> pairs_per_deme;
};
#endif
