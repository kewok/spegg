#include <thrust/gather.h>

//#include "species_specific_sampling_input.h"
#include <math/random_variables_functions.h>
#include <util/which_function.h>
#include <species/add_kids/sampling_input_mating.h>

/* Implementation of underlying methods for the different species */

#define MAGIC_THRESHOLD 50 // specifies when to switch from poisson to gaussian
#define MALE_CODE 1

void SamplingInput_Mating::determine_mate_sampling_scheme(int species_ID)
	{
	std::cout << "Try to find the sampling scheme" << std::endl;
	const char *filename = "deme_config.txt";

	Config cfg;
	
	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No deme config file." << std::endl;
		exit(1);
		}
	catch(const ParseException &pex)
		{
		std::cerr << "Your " << pex.getFile() << " file is incorrectly specified. Make sure you check on or about line: " << pex.getLine() << " - " << pex.getError() << std::endl;
		exit(1);
		}

	const Setting& root = cfg.getRoot();

	// Read in the input

	const Setting &species_specification = root["species_data"];

	bool does_sampling_scheme_exist = species_specification[species_ID].lookupValue("MATE_SAMPLING_SCHEME", sampling_scheme);

	std::cout << "sampling scheme is: " << sampling_scheme << std::endl;

	if (!does_sampling_scheme_exist)
		{
		std::cerr << "No 'MATE_SAMPLING_SCHEME' value in species data in deme_config.txt. Simulation will exit." << std::endl;
		exit(1);
		}
	}


SamplingInput_Mating::SamplingInput_Mating(Assortative_mating_parents *mating_parents, int Sampling_Parent)
	{
	/* If females are doing the sampling */
	if (Sampling_Parent != MALE_CODE)
		{
		list_of_individuals_potentially_conducting_sampling.resize(mating_parents->female_parents.size());
		list_of_individuals_potentially_subject_to_sampling.resize(mating_parents->male_parents.size());

		list_of_individuals_potentially_conducting_sampling = mating_parents->female_parents;
		list_of_individuals_potentially_subject_to_sampling = mating_parents->male_parents;
		sampleable_individuals_per_deme.resize(mating_parents->Num_Demes);
		sampleable_individuals_per_deme = mating_parents->reproductive_males_per_deme;
		}

	else
		{
		list_of_individuals_potentially_conducting_sampling.resize(mating_parents->male_parents.size());
		list_of_individuals_potentially_subject_to_sampling.resize(mating_parents->female_parents.size());

		list_of_individuals_potentially_conducting_sampling = mating_parents->male_parents;
		list_of_individuals_potentially_subject_to_sampling = mating_parents->female_parents;
		sampleable_individuals_per_deme.resize(mating_parents->Num_Demes);
		sampleable_individuals_per_deme = mating_parents->reproductive_females_per_deme;
		}

	determine_mate_sampling_scheme(mating_parents->species_ID);

	deme_affiliation_of_sampling_individuals.resize(list_of_individuals_potentially_conducting_sampling.size());
	
	thrust::gather(list_of_individuals_potentially_conducting_sampling.begin(), 
		       list_of_individuals_potentially_conducting_sampling.end(), 
		       mating_parents->deme.begin(),
		       deme_affiliation_of_sampling_individuals.begin());

	mating_scheme = mating_parents->demeParameters->species_specific_values["mating_scheme"];
	
	Num_Demes = mating_parents->Num_Demes;
	}

void SamplingInput_Mating::determine_number_of_individuals_sampled(Assortative_mating_parents *mating_parents)
	{
	if (mating_scheme == 0.0)
		determine_number_of_individuals_to_be_sampled_poisson( mating_parents );

	if (mating_scheme == 1.0)
		determine_number_of_individuals_to_be_sampled_fixed( mating_parents ) ;
	}


// Methods to determine how many individuals each individual samples
void SamplingInput_Mating::determine_number_of_individuals_to_be_sampled_poisson(Assortative_mating_parents *mating_parents) 
	{
	int number_of_individuals_doing_the_sampling = list_of_individuals_potentially_conducting_sampling.size();

	number_of_other_individuals_sampled.resize(number_of_individuals_doing_the_sampling);

	// For each individual, based on their deme, identify the expected number of other individuals they are likely to sample
	thrust::device_vector<float> mean_numbers_of_others_sampled( number_of_individuals_doing_the_sampling );
	
	thrust::gather(deme_affiliation_of_sampling_individuals.begin(), 
		       deme_affiliation_of_sampling_individuals.end(), 
		       mating_parents->demeParameters->get_vector_ptr("mean_number_of_others_sampled"),
		       mean_numbers_of_others_sampled.begin());

	// draw a poisson, where mean_numbers_of_others_sampled is the mean.
	draw_poisson_different_parameters(number_of_individuals_doing_the_sampling, 
					  mean_numbers_of_others_sampled, 
					  number_of_other_individuals_sampled,
					  mating_parents->gen);

	// Make sure you sample no more than there are individuals in the other species deme by:
		// 1. Apply correction for small deme size: cannot sample more individuals than there are individuals in the deme.
	thrust::device_vector<int> maximum_sampleable_individuals( number_of_individuals_doing_the_sampling  );
	thrust::gather(deme_affiliation_of_sampling_individuals.begin(),
		       deme_affiliation_of_sampling_individuals.end(),
		       sampleable_individuals_per_deme.begin(),
		       maximum_sampleable_individuals.begin());

	thrust::transform( number_of_other_individuals_sampled.begin(),  number_of_other_individuals_sampled.end(),  maximum_sampleable_individuals.begin(),  number_of_other_individuals_sampled.begin(), thrust::minimum<int>());
	}


// Method to set the same number of others sampled for all individuals in a deme
void SamplingInput_Mating::determine_number_of_individuals_to_be_sampled_fixed(Assortative_mating_parents *mating_parents) 
	{
	int number_of_individuals_doing_the_sampling = list_of_individuals_potentially_conducting_sampling.size();

	number_of_other_individuals_sampled.resize(number_of_individuals_doing_the_sampling);

	// Assume all individuals sample the same number of individuals in the other species.

	// specify the number of other individuals sampled, which we assume is set deterministically
	thrust::gather(deme_affiliation_of_sampling_individuals.begin(), 
		       deme_affiliation_of_sampling_individuals.end(), 
		       mating_parents->demeParameters->get_vector_ptr("mean_number_of_others_sampled"),
		       number_of_other_individuals_sampled.begin());
	

	// Apply correction for small deme size: cannot sample more individuals than there are individuals in the deme.
	thrust::device_vector<int> maximum_sampleable_individuals( number_of_individuals_doing_the_sampling  );

	thrust::gather(deme_affiliation_of_sampling_individuals.begin(),
		       deme_affiliation_of_sampling_individuals.end(),
		       sampleable_individuals_per_deme.begin(),
		       maximum_sampleable_individuals.begin());

	thrust::transform( number_of_other_individuals_sampled.begin(),  number_of_other_individuals_sampled.end(),  maximum_sampleable_individuals.begin(),  number_of_other_individuals_sampled.begin(), thrust::minimum<int>());
	}
//////////////////////
///
/// Repeat for species 1, 2, 3, 4, etc...
///
/////////////////////
