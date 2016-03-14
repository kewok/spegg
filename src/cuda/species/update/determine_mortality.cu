#include <species/update/updatebehavior.h>


void UpdateBehavior::determine_mortality(inds_stochastic *species)
	{
	/*
	* Code that actually changes the individual's vital state (dead=0, alive=1) during each time step. The value of the individual's probability of survival (stored in \code species->phenotype[MORTALITY_PHENOTYPE] \endcode is set in an external function (e.g., something like \code update_mySpecies::determine_survivorship_probability() \endcode or something like that.).
	*/

	//Specify the individuals indices
	thrust::host_vector<int> individuals(species->size);
	thrust::sequence(individuals.begin(), individuals.begin() + species->size, 0);

	// Draw the random numbers
	thrust::host_vector<float> rand(species->size);
	for (int i=0; i < rand.size(); i++)
		{
		rand[i] = gsl_rng_uniform(species->gen);
		}
	float *rand_ptr = (&rand[0]);

	//Set up mortality.
	simulate_mortality mortality_functor(rand_ptr);

	//Perform mortality operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals.begin(), species->status.begin(), species->phenotype[species->MORTALITY_PHENOTYPE_INDEX].begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + species->size, species->status.begin() + species->size, species->phenotype[species->MORTALITY_PHENOTYPE_INDEX].begin() + species->size)),
			 mortality_functor);
	}
