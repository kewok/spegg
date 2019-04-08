#ifndef SAMPLING_INPUT_MATING_H
#define SAMPLING_INPUT_MATING_H

#include <iostream>
#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <util/sampling_input.h>
#include <species/add_kids/assortative_mating_parents_class.h>

class SamplingInput_Mating : public SamplingInput
	{
	/*
	* A class whereby females or males of a given species will "sample" individuals of the opposite sex and evaluate their suitability as mates. A common example for a situation where this class would be necessary would be when female preference varies across individual females, who must evaluate a sequence of males and select a mate according to their innate preference function.
	*/
	public:
		SamplingInput_Mating(Assortative_mating_parents *mating_parents, int Sampling_Parent);

		void determine_number_of_individuals_sampled(Assortative_mating_parents *mating_parents);

		float mating_scheme;
	protected:
		// Stock algorithms for determining the number of individuals sampled
		void determine_number_of_individuals_to_be_sampled_fixed(Assortative_mating_parents *mating_parents); 
		void determine_number_of_individuals_to_be_sampled_poisson(Assortative_mating_parents *mating_parents); 

		void determine_mate_sampling_scheme(int species_index);
	};
#endif
