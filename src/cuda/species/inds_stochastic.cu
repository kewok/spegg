#include <species/inds_stochastic.h>

inds_stochastic::inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds(size_val, maxsize_val, ndemes, species_ID_val)
	{
/*
*
* Initialize the curand generator gen based on the seed argument, using CURAND_RNG_PSEUDO_DEFAULT. Draw 100 random uniform variables and store them in the vector rand; then discard the contents of the rand vector. This is all done to initialize the CUDA random number generator object. Note the rand vector will then be deallocated once the prime_random_number_generator scope is ended. For reasons not entirely clear, this can't seem to be done inside an external function prime_random_number_generator
*
*/
	seed = seed_val;

	int size = 100;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, seed);

	//curand declarations
	thrust::device_vector<float> rand(size);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(gen, rand_ptr, size); // priming up the random number generator takes some time, get it done early.
	rand.clear();

	//Specify the indices among the phenotypes for the fitness components
	MORTALITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
	FECUNDITY_PHENOTYPE_INDEX = (int) demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	}

