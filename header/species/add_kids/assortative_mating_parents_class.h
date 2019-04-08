#ifndef ASSOCIATIVE_MATING_PARENTS_CLASS_H
#define ASSOCIATIVE_MATING_PARENTS_CLASS_H

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <curand.h>
#include <species/inds.h>
#include <species/add_kids/parents_class.h>
#include <math.h>
class Assortative_mating_parents : public Parents
	{
	public:
		Assortative_mating_parents(inds_stochastic *species);

		void determine_parent_pair_probability(thrust::device_vector<float> *&phenotypes);
		thrust::device_vector<float> parental_pair_probability;

		void draw_parents();

		/* female_parents, male_parents differs from females_list, males_list. In the assortative mating context, female_parents, male_parents merely stores the indices of individuals that can potentially be parents, while females_list,males_list stores the potential sampling space (ie., has duplicates) */
		thrust::device_vector<int> females_list; 
		thrust::device_vector<int> males_list;

		int total_reproductive_females;
		int total_reproductive_males;
	
		thrust::device_vector<int> pair_demes;

		int species_ID;
	protected:

		int ASSORTATIVE_MATING_PHENOTYPE_INDEX;

		int mate_sampling_scheme;

		thrust::device_vector<int> number_of_males_sampled;
		thrust::device_vector<float> assortative_mating_trait;

		void Setup_ParentPairProbabilities();

		void Generate_Parents_List();

		void Generate_Females_List();
		void Generate_Males_List();
	};


struct assign_males
	{
	int *cumulative_males, *number_of_males_by_deme;
	assign_males(int *cumul_males, int *num_males_in_deme) : cumulative_males (cumul_males), number_of_males_by_deme(num_males_in_deme)
	{};
	
	/*
		Elements in the tuple.
		----------------------
		0: index
		1: female population
		2: index of male
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int females_pop = thrust::get<1>(t);
		thrust::get<2>(t) = (thrust::get<0>(t) - cumulative_males[females_pop])%(number_of_males_by_deme[females_pop]) + cumulative_males[females_pop];
		}
	};

struct pairwise_mating_probability
	{
	float *assortative_mating_trait;
	float *assortative_mating_value;
	pairwise_mating_probability(float *assort_mating_trait, float *assort_mating_value) : assortative_mating_trait (assort_mating_trait), assortative_mating_value (assort_mating_value)
	{};

	/*
		Elements in the tuple.
		----------------------
		0: female index
		1: male index
		2: deme the parents come from
		3: female reproductive potential
		4: pairwise mating probability
	*/
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		float female_assortative_mating_trait;
		float male_assortative_mating_trait;
	
		female_assortative_mating_trait = assortative_mating_trait[thrust::get<0>(t)];
		male_assortative_mating_trait = assortative_mating_trait[thrust::get<1>(t)];

		float importance_of_assortative_mating = assortative_mating_value[thrust::get<2>(t)];

		/* the relative probability that male thrust::get<1>(t) mates with female thrust::get<0>(t) is just the inverse of the difference between their assortative mating traits (i.e., the further these values are the less likely they are to mate) */

		thrust::get<4>(t) = exp(-pow(female_assortative_mating_trait - male_assortative_mating_trait,2)/ importance_of_assortative_mating) * thrust::get<3>(t);
		}
	};


#endif
