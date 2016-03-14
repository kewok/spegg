#ifndef PARENTS_CLASS_H
#define PARENTS_CLASS_H

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <curand.h>
#include <species/inds_stochastic.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <util/thrust_functors.h>


class Parents
	{
	friend class SamplingInput_Mating;

	public:
		Parents(inds_stochastic *species);
	
		/* since the arguments passed to these functions should NOT change, pure virtual is OK to use here: */
		virtual void setup_parents();

		thrust::device_vector<int> female_parents;
		thrust::device_vector<int> male_parents;

		int Potential_Number_of_Kids;
		thrust::device_vector<int> reproductive_potential_per_deme;

		thrust::device_vector<int> reproductive_males_per_deme;
		thrust::device_vector<int> reproductive_females_per_deme;

		thrust::device_vector<int> kids_per_mom;

		/* stuff related to selection as parents */
		thrust::device_vector<float> probability_individual_becomes_female_parent;
		thrust::device_vector<float> probability_individual_becomes_male_parent;

	protected:
		curandGenerator_t gen;

		DemeSettings *demeParameters;

		int size;
		int Num_Demes;
		int Total_Number_of_Parents;

		thrust::device_vector<float> *phenotype;
		thrust::device_vector<int> deme;
		thrust::device_vector<int> sex;

		/* stuff related to eligibility */
		thrust::device_vector<int> will_reproduceF;
		thrust::device_vector<int> will_reproduceM;

		virtual void finalize_parental_reproductive_probabilities();
		void determine_parental_reproductive_potential();
		virtual void determine_female_parent_eligibility();
		virtual void determine_male_parent_eligibility();
		virtual void determine_probability_individual_becomes_female_parent();
		virtual void determine_probability_individual_becomes_male_parent();
		virtual void female_fecundity();
	};

#endif
