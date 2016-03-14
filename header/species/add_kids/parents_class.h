#ifndef PARENTS_CLASS_H
#define PARENTS_CLASS_H

#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <species/inds_stochastic.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <util/thrust_functors.h>

class Parents
	{
	friend class SamplingInput_Mating;

	public:
		Parents(inds_stochastic *species);

		virtual void setup_parents();

		thrust::host_vector<int> female_parents;
		thrust::host_vector<int> male_parents;

		int Potential_Number_of_Kids;
		thrust::host_vector<int> reproductive_potential_per_deme;

		thrust::host_vector<int> reproductive_males_per_deme;
		thrust::host_vector<int> reproductive_females_per_deme;

		thrust::host_vector<int> kids_per_mom;

		/* stuff related to selection as parents */
		thrust::host_vector<float> probability_individual_becomes_female_parent;
		thrust::host_vector<float> probability_individual_becomes_male_parent;

	protected:
		gsl_rng * gen;

		DemeSettings *demeParameters;

		int size;
		int Num_Demes;
		int Total_Number_of_Parents;
	
		thrust::host_vector<float> *phenotype;
		thrust::host_vector<int> deme;
		thrust::host_vector<int> sex;
	
		/* stuff related to eligibility */
		thrust::host_vector<int> will_reproduceF;
		thrust::host_vector<int> will_reproduceM;

		virtual void finalize_parental_reproductive_probabilities();
		void determine_parental_reproductive_potential();
		virtual void determine_female_parent_eligibility();
		virtual void determine_male_parent_eligibility();
		virtual void determine_probability_individual_becomes_female_parent();
		virtual void determine_probability_individual_becomes_male_parent();
		virtual void female_fecundity();
	};

#endif
