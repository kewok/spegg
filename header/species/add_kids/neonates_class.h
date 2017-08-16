#ifndef EGGS_NEONATES_H
#define EGGS_NEONATES_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <species/inds_stochastic.h>
#include <species/deme_specific_data_class.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <util/amplify.h>
#include <math/mating_thrust_prob_table_demes.h>
#include <math/random_variables_functions.h>

class EggsNeonates
	{
	public:
		EggsNeonates(inds_stochastic *species, thrust::host_vector<int> &kids_per_mom);

		void inherit_genotypes(thrust::host_vector<float> &probability_individuals_become_mothers,
					 thrust::host_vector<float> &probability_individuals_become_fathers);

		int previous_pop_size;
		int Total_Number_of_Neonates;

		thrust::host_vector<int> Neonates_per_Deme;

		thrust::host_vector<int> kids_deme; 

		protected:
			class inds *species;
			gsl_rng *gen;	
			int nloci, nphen, Num_Demes;

			thrust::host_vector<float> mutation_magnitude;
			thrust::host_vector<float> mutation_rate;

			void Determine_Neonate_Population_Sizes(DemeSettings *subpopParameters,
								thrust::host_vector<int> &everybodys_deme,
								thrust::host_vector<int> &kids_per_mom,			
								thrust::host_vector<int> &current_deme_sizes,
								thrust::host_vector<int> &maximum_deme_sizes);

			void get_maternally_derived_genotype(thrust::host_vector<float> &probability_individuals_become_mothers, thrust::host_vector<float> *&mgenotype, thrust::host_vector<float> *&fgenotype);

			void get_paternally_derived_genotype(thrust::host_vector<float> &probability_individuals_become_fathers, thrust::host_vector<float> *&mgenotype, thrust::host_vector<float> *&fgenotype);

			void egg_mortality(DemeSettings *subpopParameters);

			void recombine(thrust::host_vector<float> &rand,
						 thrust::host_vector<int> &parent,
						 thrust::host_vector<int> &parity,
						 thrust::host_vector<float> *&parents_fgenotype,
						 thrust::host_vector<float> *&parents_mgenotype,
						 thrust::host_vector<float> &kids_genotype,
						 int locus_ID);

			thrust::host_vector<float> recomb_rate;

			void mutate(thrust::host_vector<float> *&parents_fgenotype, thrust::host_vector<float> *&parents_mgenotype);

			void prepare_genotype_phenotype_map();

			void integrate_kids();
	};

// make sure there are no more kids than spaces available
struct adjust_kids_functor
	{
	/*
		Elements in the tuple.
		----------------------
		0: subpopulation kid size
		1: subpopulation size already there
		2: subpopulation carrying capacity
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
			if (thrust::get<0>(t) + thrust::get<1>(t) > thrust::get<2>(t)) 
			{
			if (thrust::get<2>(t) - thrust::get<1>(t) >= 0 )
					{
					thrust::get<0>(t) = thrust::get<2>(t) - thrust::get<1>(t);
					}
			if (thrust::get<2>(t) - thrust::get<1>(t) < 0 ) // if too crowded, no kids.
					{
					thrust::get<0>(t) = 0;
					}
			}
		}
	};

struct recombination_functor
	{
	float *fgenotype, *mgenotype;
	float recomb_rate;
	recombination_functor(float *fgene, float *mgene, float rate) : fgenotype(fgene), mgenotype(mgene), recomb_rate(rate)
	{};
	
	/*
		Elements in the tuple.
		----------------------
		0: parent index
		1: random number
		2: parity
		3: genotype
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int parent_index = thrust::get<0>(t);

		if (thrust::get<1>(t) < recomb_rate) 
			{
			thrust::get<2>(t) = (thrust::get<2>(t) + 1)%2;
			}
		if (thrust::get<2>(t) == 0)
			thrust::get<3>(t) = fgenotype[parent_index];
		if (thrust::get<2>(t) == 1)
			thrust::get<3>(t) = mgenotype[parent_index];

		float answer = thrust::get<3>(t);
		}
	};

#endif
