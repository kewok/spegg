#include <species/inds_stochastic.h>

inds_stochastic::inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds(size_val, maxsize_val, ndemes, species_ID_val)
	{
/*
*
* Initialize the random number generator gen based on the seed argument
*
*/
	seed = seed_val;

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	gen = gsl_rng_alloc (T);
	gsl_rng_set(gen, seed);
	
	//Specify the indices among the phenotypes for the fitness components
	MORTALITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
	FECUNDITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	}

