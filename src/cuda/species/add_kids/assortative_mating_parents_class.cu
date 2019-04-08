#include <species/inds.h>
#include <species/deme_specific_data_class.h>
#include <species/add_kids/parents_class.h>
#include <species/add_kids/sampling_input_mating.h>
#include <species/add_kids/assortative_mating_parents_class.h>

#include <util/reduce_by_key_with_zeroes.h>
#include <util/thrust_functors.h>
#include <util/amplify.h>
#include <util/which_function.h>
#include <util/gather_values_by_deme.h>
#include <util/footimer2.h>
#include <util/sampling_event.h>
#include <math/mating_thrust_prob_table_demes.h>


#define SAMPLING_PARENT_SEX 0 // Assume females do the sampling

Assortative_mating_parents::Assortative_mating_parents(inds_stochastic *species) : Parents(species)
	{
	determine_female_parent_eligibility();
	determine_male_parent_eligibility();

	female_parents.resize(size);
	male_parents.resize(size);	

	thrust::fill(female_parents.begin(), female_parents.end(), -1);
	thrust::fill(male_parents.begin(), male_parents.end(), -1);

	reduce_by_key_with_zeros(deme, will_reproduceM, reproductive_males_per_deme, size, Num_Demes); 
	reduce_by_key_with_zeros(deme, will_reproduceF, reproductive_females_per_deme, size, Num_Demes);
	
	/* determine who the potential male, female parents are */
	total_reproductive_females = thrust::reduce(reproductive_females_per_deme.begin(), reproductive_females_per_deme.end());
	total_reproductive_males = thrust::reduce(reproductive_males_per_deme.begin(), reproductive_males_per_deme.end());

	female_parents.resize(total_reproductive_females);
	male_parents.resize(total_reproductive_males);

	which_equal_to(will_reproduceF, female_parents, 1);
	which_equal_to(will_reproduceM, male_parents, 1);

	// TODO: Add an error handler that will report an error if assortative_mating_phenotyope_index is absent at this stage
	ASSORTATIVE_MATING_PHENOTYPE_INDEX = (int) species->demeParameters->species_specific_values["ASSORTATIVE_MATING_PHENOTYPE"];

	assortative_mating_trait.resize(size);
	thrust::copy(species->phenotype[ASSORTATIVE_MATING_PHENOTYPE_INDEX].begin(), species->phenotype[ASSORTATIVE_MATING_PHENOTYPE_INDEX].begin() + size, assortative_mating_trait.begin());

	species_ID = species->species_ID;
	mate_sampling_scheme = 0;
	}

void Assortative_mating_parents::determine_parent_pair_probability(thrust::device_vector<float> *&phenotypes)
	{

	Generate_Parents_List();
/*
	Generate_Females_List();
	Generate_Males_List();
*/
	// Get the deme in which the pair is found
	pair_demes.resize(females_list.size());
	thrust::gather(females_list.begin(), females_list.end(), deme.begin(), pair_demes.begin());
	parental_pair_probability.resize(females_list.size());

	float *assortative_mating_trait_ptr = raw_pointer_cast(&assortative_mating_trait[0]);
	float *assortative_mating_value_ptr = raw_pointer_cast(demeParameters->get_vector_ptr("assortative_mating_value"));
	
	pairwise_mating_probability probability_pairs_mate(assortative_mating_trait_ptr, assortative_mating_value_ptr);

	thrust::device_vector<float> female_sizes(females_list.size());
	thrust::gather(females_list.begin(), females_list.end(), phenotypes[1].begin(), female_sizes.begin());

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(females_list.begin(), males_list.begin(), pair_demes.begin(), female_sizes.begin(), parental_pair_probability.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(females_list.end(), males_list.end(), pair_demes.end(), female_sizes.end(), parental_pair_probability.end())),
			 probability_pairs_mate);
	}


void Assortative_mating_parents::Generate_Parents_List()
	{
	SamplingInput_Mating *PotentialParents;
	PotentialParents = new SamplingInput_Mating(this, SAMPLING_PARENT_SEX);
	PotentialParents->determine_number_of_individuals_sampled(this);

	SamplingEvent *sampling_event;
	sampling_event = sampling_event->create_SamplingEvent(PotentialParents, gen);
	sampling_event->sample();

	females_list.resize(sampling_event->sampling_individuals_indices.size());
	males_list.resize(sampling_event->sampling_individuals_indices.size());

	thrust::copy(sampling_event->sampling_individuals_indices.begin(), sampling_event->sampling_individuals_indices.end(), females_list.begin());

	thrust::copy(sampling_event->sampled_individuals_indices.begin(), sampling_event->sampled_individuals_indices.end(), males_list.begin());

	delete PotentialParents;
	delete sampling_event;
	}

/******* 
*
*
These routines are for the brute-force, calculate all pairwise interaction routines 
*
*
*******/

void Assortative_mating_parents::Generate_Females_List()
	{
	/* flay out the list of females, in essence, this is the R operation rep(female_parents, each=reproductive_males_per_deme) */

	thrust::device_vector<int> females_deme(total_reproductive_females);
	thrust::gather(female_parents.begin(), female_parents.begin() + total_reproductive_females, deme.begin(), females_deme.begin());

	/* determine how many times each female repreats */
	thrust::device_vector<int> females_repeat(total_reproductive_females);
	thrust::gather(females_deme.begin(), females_deme.end(), reproductive_males_per_deme.begin(), females_repeat.begin());

	amplify(female_parents, females_repeat,  females_list);
	}

void Assortative_mating_parents::Generate_Males_List()
	{
	/* now do males. In essence, this is the R operation rep(male_parents, reproductive_females_per_deme) */
	thrust::device_vector<int> females_deme(males_list.size());
	thrust::gather(females_list.begin(), females_list.end(), deme.begin(), females_deme.begin());

	thrust::device_vector<int> cumulative_males_by_deme(Num_Demes);
	thrust::exclusive_scan(reproductive_males_per_deme.begin(), reproductive_males_per_deme.end(), cumulative_males_by_deme.begin());
	
	thrust::device_vector<int> indices(females_list.size());
	thrust::sequence(indices.begin(), indices.end(), 0);
	
	males_list.resize(females_list.size());
	
	int *cumulative_males_by_deme_ptr = raw_pointer_cast(&cumulative_males_by_deme[0]);
	int *males_per_deme_ptr = raw_pointer_cast(&reproductive_males_per_deme[0]);
	
	assign_males male_vals(cumulative_males_by_deme_ptr, males_per_deme_ptr);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(indices.begin(), females_deme.begin(), males_list.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(indices.end(), females_deme.end(), males_list.end())),
			 male_vals);
	
	thrust::gather(males_list.begin(), males_list.end(), male_parents.begin(), males_list.begin());
	}

