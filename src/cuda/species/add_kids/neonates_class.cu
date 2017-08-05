#include <species/add_kids/neonates_class.h>

/* One problem with a separate neonates class is that all this stuff has to get copied back into example species. Performance cost seems minimal but more extensive profiling seems warranted. */

EggsNeonates::EggsNeonates(inds_stochastic *species, thrust::device_vector<int> &kids_per_mom) 
	{
	this->species = species;
	this->gen = species->gen;

	nphen = species->nphen;
	nloci = species->nloci;

	//Copy in the recombination map
	recomb_rate.resize(nloci);
	thrust::copy(species->demeParameters->GeneticArchitecture->recombination_rates.begin(), species->demeParameters->GeneticArchitecture->recombination_rates.begin() + nloci, recomb_rate.begin());

	previous_pop_size = species->size;

	Num_Demes = species->Num_Demes;
	Neonates_per_Deme.resize(Num_Demes);
	thrust::fill(Neonates_per_Deme.begin(), Neonates_per_Deme.end(), 0);

	Determine_Neonate_Population_Sizes(species->demeParameters, species->deme, kids_per_mom, species->deme_sizes, species->max_deme_sizes);

	kids_deme.resize(Total_Number_of_Neonates);

	amplify_sequence( Neonates_per_Deme, Num_Demes, kids_deme );

	// See to it that the allocated deme gets the neonates assigned.
	integrate_kids();

	mutation_magnitude.resize(Total_Number_of_Neonates);
	mutation_rate.resize(Total_Number_of_Neonates);
	}

void EggsNeonates::Determine_Neonate_Population_Sizes(DemeSettings *demeParameters,
						      thrust::device_vector<int> &everybodys_deme,
						      thrust::device_vector<int> &kids_per_mom,
						      thrust::device_vector<int> &current_deme_sizes,
						      thrust::device_vector<int> &maximum_deme_sizes)
	{
	/*
	* Because there are often more neonates than max_deme_sizes, this function culls the surplus neonates at random by determining the number of neonates each deme can contribute
	*/
	cudaThreadSynchronize();
	reduce_by_key_with_zeros(everybodys_deme, kids_per_mom, Neonates_per_Deme, previous_pop_size, Num_Demes); 

	// Make sure no population has more kids than there are spaces
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(Neonates_per_Deme.begin(), current_deme_sizes.begin(), maximum_deme_sizes.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(Neonates_per_Deme.end(), current_deme_sizes.end(), maximum_deme_sizes.end())),
	                 adjust_kids_functor());

	Total_Number_of_Neonates = thrust::reduce(Neonates_per_Deme.begin(), Neonates_per_Deme.end());
	}


void EggsNeonates::inherit_genotypes(thrust::device_vector<float> &probability_individuals_become_mothers,
				 thrust::device_vector<float> &probability_individuals_become_fathers)
	{
	get_maternally_derived_genotype(probability_individuals_become_mothers, species->mgenotype, species->fgenotype);
	get_paternally_derived_genotype(probability_individuals_become_fathers, species->mgenotype, species->fgenotype);
	mutate(species->fgenotype, species->mgenotype);
	}


void EggsNeonates::get_maternally_derived_genotype(thrust::device_vector<float> &probability_individuals_become_mothers,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype)
	{
	/*
	* Calculate the offspring's genotypes at the maternally inherited loci. This function will determine who the mother is, and generate a haploid gamete from the mother that will be grafted on to the offspring.
	*/
	mating_ThrustProbTable_demes at;
	thrust::device_vector<int> mother_index(Total_Number_of_Neonates);
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
/*
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine mothers.
*/
	at.setup(probability_individuals_become_mothers.begin(), probability_individuals_become_mothers.begin() + previous_pop_size);
	curandGenerateUniform(gen, rand_ptr, Total_Number_of_Neonates);
	at.determine_key_offsets( Num_Demes, species->deme_sizes );
	at.adjust_randoms(rand.begin(), rand.end(), kids_deme.begin(), kids_deme.end());

	at.draw(rand.begin(), rand.end(), mother_index.begin());
	thrust::copy(mother_index.begin(), mother_index.end(), species->maternal_id.begin() + previous_pop_size);

//Initialize parity to zeroes
//Parity vector is used to keep track of where the recombination is happening.
	thrust::fill(parity.begin(), parity.end(), 0);


//Recombination for fgenotype
	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(gen, rand_ptr, Total_Number_of_Neonates);
        	recombine(rand, mother_index, parity, fgenotype, mgenotype, fgenotype[i], i);
		}
	}

void EggsNeonates::get_paternally_derived_genotype(thrust::device_vector<float> &probability_individuals_become_fathers,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype)
	{
	/*
	* Calculate the offspring's genotypes at the paternally inherited loci. This function will determine who the mother is, and generate a haploid gamete from the mother that will be grafted on to the offspring.
	*/
	mating_ThrustProbTable_demes at;
	thrust::device_vector<int> father_index(Total_Number_of_Neonates);
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	
	/*
	Directly use phenotype[1] as the reproductive probability.
	Reproductive probability = phenotype[1].
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine fathers.
	*/
	
	at.setup(probability_individuals_become_fathers.begin(), probability_individuals_become_fathers.begin() + previous_pop_size);
	curandGenerateUniform(gen, rand_ptr, Total_Number_of_Neonates);

	at.determine_key_offsets( Num_Demes, species->deme_sizes );	
	at.adjust_randoms(rand.begin(), rand.end(), kids_deme.begin(), kids_deme.end());

	at.draw(rand.begin(), rand.end(),father_index.begin());
	thrust::copy(father_index.begin(), father_index.end(), species->paternal_id.begin() + previous_pop_size);

	//Reset parity to zeroes
	thrust::fill(parity.begin(), parity.end(), 0);

	//Recombination for mgenotype
	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(gen, rand_ptr, Total_Number_of_Neonates);
		recombine(rand, father_index, parity, fgenotype, mgenotype, mgenotype[i], i);
		}
	}

void EggsNeonates::mutate(thrust::device_vector<float> *&mgenotype, thrust::device_vector<float> *&fgenotype)
	{
	/*
	* Determine, for each offspring's locus, whether there will be a mutation at that locus, and if so what the magnitude of that mutation will be and how that changes the offspring's allelic value. The current behavior assumes mutations are gaussian about the parental allelic value, and that the mutation parameters (mutation rate and sd of gaussian) vary by deme. Possible expansions include allowing the mutation rate itself to be an individual-specific phenotype, or alternative mutational models (e.g., point mutations that have a categorical rather than quantitative effect.
	*/
	thrust::identity<int> identity;

	// Note throughout we use 2*Total_Number_of_Neonates; this is for performance reasons. drawing gaussians can be expensive so we want to parallelize as much of this as we can.
	thrust::device_vector<float> mutation_size(2*Total_Number_of_Neonates); 

	thrust::device_vector<float> mutation_prob(2*Total_Number_of_Neonates);

	thrust::device_vector<int> mutation_event(2*Total_Number_of_Neonates);
	thrust::fill(mutation_event.begin(), mutation_event.begin() + 2*Total_Number_of_Neonates, 0);
	// no mutation at the original sex determining locus

	for (int i = 0 ; i < nloci ; i++) 
		{
		// assign, for each offspring, the expected magnitude of the mutation (sd) and the expected mutation rate according to the their deme. At some point the performance of this regime needs to be profiled, because there are a lot of steps involved here.

		thrust::gather(kids_deme.begin(), kids_deme.begin() + Total_Number_of_Neonates, species->demeParameters->GeneticArchitecture->get_mutation_magnitudes_ptr(i), mutation_magnitude.begin());

		thrust::device_vector<float> mutation_size(2*Total_Number_of_Neonates); 

		thrust::copy(mutation_magnitude.begin(), mutation_magnitude.begin() + Total_Number_of_Neonates, mutation_size.begin());
		thrust::copy(mutation_magnitude.begin(), mutation_magnitude.begin() + Total_Number_of_Neonates, mutation_size.begin() + Total_Number_of_Neonates);

		thrust::gather(kids_deme.begin(), kids_deme.begin() + Total_Number_of_Neonates, species->demeParameters->GeneticArchitecture->get_mutation_rates_ptr(i), mutation_rate.begin());

		thrust::copy(mutation_rate.begin(), mutation_rate.begin() + Total_Number_of_Neonates, mutation_prob.begin());
		thrust::copy(mutation_rate.begin(), mutation_rate.begin() + Total_Number_of_Neonates, mutation_prob.begin() + Total_Number_of_Neonates);

		thrust::device_vector<float> mutation(2*Total_Number_of_Neonates);

		// Determine whether mutations will occur
		draw_bernoulli_different_parameters(2*Total_Number_of_Neonates, mutation_prob, mutation_event, gen);
		// Assume the mutation is symmetric about the current allelic value:
		draw_gaussian_different_parameters(2*Total_Number_of_Neonates, (float) 0.0, mutation_size, mutation, gen);
	
		// If there is a mutation, add it to the offspring's 
		thrust::transform_if(mgenotype[i].begin() + previous_pop_size, mgenotype[i].begin()  + previous_pop_size + Total_Number_of_Neonates, mutation.begin(), mutation_event.begin(), mgenotype[i].begin() + previous_pop_size, thrust::plus<float>(), identity);
		thrust::transform_if(fgenotype[i].begin() + previous_pop_size, fgenotype[i].begin()  + previous_pop_size + Total_Number_of_Neonates, mutation.begin() + Total_Number_of_Neonates, mutation_event.begin() + Total_Number_of_Neonates, fgenotype[i].begin()  + previous_pop_size, thrust::plus<float>(), identity);
		}
	}


void EggsNeonates::recombine(thrust::device_vector<float> &rand,
			     thrust::device_vector<int> &parent,
			     thrust::device_vector<int> &parity,
			     thrust::device_vector<float> *&parents_fgenotype,
			     thrust::device_vector<float> *&parents_mgenotype,
			     thrust::device_vector<float> &kids_genotype,
			     int locus_ID)
	{
/*
* A method that implements recombination. Note that in many cases, the arguments for either (parents_fgenotype or parents_mgenotype) and kids_genotype will be the same vector. This function increments the kids_genotype to begin at the (previous_pop_size)th index. This should probably be refactored so that there will be a kids_mgenotype and kids_fgenotype vector that is local to the neonates class, and that later gets copied into the inds class, i.e. something like: 
\code{.cpp}
	thrust::copy(kids_fgenotype[i].begin(), kids_fgenotype[i].end(), fgenotype[i].begin() + size);
\endcode
*/
	//Set up recombination functor.
	float *fgenotype_ptr = raw_pointer_cast(&parents_fgenotype[locus_ID][0]);
	float *mgenotype_ptr = raw_pointer_cast(&parents_mgenotype[locus_ID][0]);
	recombination_functor rfunc(fgenotype_ptr, mgenotype_ptr, recomb_rate[locus_ID]);
	
	//Perform recombination with arbitrary transform.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(parent.begin(), rand.begin(), parity.begin(), kids_genotype.begin() + previous_pop_size)),
			 thrust::make_zip_iterator(thrust::make_tuple(parent.begin() + Total_Number_of_Neonates, rand.begin() + Total_Number_of_Neonates, parity.begin() + Total_Number_of_Neonates, kids_genotype.begin() +  previous_pop_size + Total_Number_of_Neonates)),
			 rfunc);
	}


void EggsNeonates::integrate_kids()
	{
	thrust::copy(kids_deme.begin(), kids_deme.begin() + Total_Number_of_Neonates, species->deme.begin() + previous_pop_size);
	thrust::sequence(species->id.begin() + previous_pop_size, species->id.begin() + previous_pop_size + Total_Number_of_Neonates, species->nextid);
	thrust::fill(species->status.begin() + previous_pop_size, species->status.begin() + previous_pop_size + Total_Number_of_Neonates, 1);
	thrust::fill(species->age.begin() + previous_pop_size, species->age.begin() + previous_pop_size + Total_Number_of_Neonates, 0);
	species->nextid += Total_Number_of_Neonates;
	species->size += Total_Number_of_Neonates;
	}
