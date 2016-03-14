#include <species/add_kids/parents_class.h>
#include <species/add_kids/mating_kernel_functors.h>
#include <util/footimer2.h>

// since we are ignoring genetics and using an external curand generator here, this is OK:

Parents::Parents(inds_stochastic *species) 
	{
	this->phenotype = species->phenotype;

	// Get relevant values from inds_stochastic object
	size = species->size;
	this->Num_Demes = species->Num_Demes;
	gen = species->gen;
	this->demeParameters = species->demeParameters;

	deme.resize(size);
	thrust::copy(species->deme.begin(), species->deme.begin() + size, deme.begin());
	
	sex.resize(size);
	thrust::copy(species->sex.begin(), species->sex.begin() + size, sex.begin());

	// Configure the parents class accordingly
	will_reproduceF.resize(size);
	will_reproduceM.resize(size);

	probability_individual_becomes_female_parent.resize(size);
	probability_individual_becomes_male_parent.resize(size);

	kids_per_mom.resize(size);
	thrust::fill(kids_per_mom.begin(), kids_per_mom.begin() + size, 0);

	reproductive_males_per_deme.resize(Num_Demes); 
	thrust::fill(reproductive_males_per_deme.begin(), reproductive_males_per_deme.end(), 0);

	reproductive_females_per_deme.resize(Num_Demes); 
	thrust::fill(reproductive_females_per_deme.begin(), reproductive_females_per_deme.end(), 0);

	reproductive_potential_per_deme.resize(Num_Demes);
	thrust::fill(reproductive_potential_per_deme.begin(), reproductive_potential_per_deme.end(), 0);
	}

void Parents::setup_parents()
	{
	determine_parental_reproductive_potential();
	finalize_parental_reproductive_probabilities();
	}

void Parents::determine_parental_reproductive_potential()
	{
	/* these two steps determine_female_parent_eligibility and determine_male_parent_eligibility are almost certainly redundant and should be refactored */
	determine_female_parent_eligibility();
	determine_male_parent_eligibility();

	reduce_by_key_with_zeros(deme, will_reproduceM, reproductive_males_per_deme, size, Num_Demes); 

	female_fecundity();

	Potential_Number_of_Kids = thrust::reduce(kids_per_mom.begin(), kids_per_mom.end());

	if (Potential_Number_of_Kids < 0)
		{
		std::cout << "Potential number of kids shouldn't be negative" << std::endl;
		for (int i=0; i < kids_per_mom.size(); i++)
			std::cout << kids_per_mom[i] << std::endl;
		}
	}

void Parents::finalize_parental_reproductive_probabilities()
	{
	determine_probability_individual_becomes_female_parent();
	determine_probability_individual_becomes_male_parent();
	}

void Parents::determine_female_parent_eligibility()
	{
	// Default template that may be overriden by derived class 
	thrust::fill(will_reproduceF.begin(), will_reproduceF.begin() + size, 0);
	
	//Set up eligibility functor.
	parental_eligibility_functor parental_eligibility(0);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), sex.begin())),
		 	thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, sex.begin() + size)),
			 parental_eligibility);	
	}

void Parents::determine_male_parent_eligibility()
	{
	// Default template that may be overriden by derived class 
	thrust::fill(will_reproduceM.begin(), will_reproduceM.begin() + size, 0);
	
	//Set up eligibility functor.
	parental_eligibility_functor parental_eligibility(1);
		
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin(), sex.begin())),
		 	thrust::make_zip_iterator(thrust::make_tuple(will_reproduceM.begin() + size, sex.begin() + size)),
			 parental_eligibility);	
	}

void Parents::female_fecundity()
	{
	int FECUNDITY_INDEX = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin(), phenotype[FECUNDITY_INDEX].begin(), probability_individual_becomes_female_parent.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(will_reproduceF.begin() + size, phenotype[FECUNDITY_INDEX].begin() + size, probability_individual_becomes_female_parent.begin() + size)),
			 female_fecundity_functor());

	thrust::copy(probability_individual_becomes_female_parent.begin(), probability_individual_becomes_female_parent.end(), kids_per_mom.begin());
	}

void Parents::determine_probability_individual_becomes_female_parent()
	{
	/* Virtual void placeholders to be replaced by species-specific parents class*/
	}

void Parents::determine_probability_individual_becomes_male_parent()
	{
	/* Virtual void placeholder to be replaced by species-specific parents class*/
	}

