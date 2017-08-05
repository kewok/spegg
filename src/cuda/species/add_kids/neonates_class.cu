#include <species/add_kids/neonates_class.h>

/* One problem with a separate neonates class is that all this stuff has to get copied back into example species. It might just not be worth it depending on the performance cost. */

EggsNeonates::EggsNeonates(inds_stochastic *species, thrust::host_vector<int> &kids_per_mom) 
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

	integrate_kids();
	mutation_magnitude.resize(Total_Number_of_Neonates);
	mutation_rate.resize(Total_Number_of_Neonates);
	}

void EggsNeonates::Determine_Neonate_Population_Sizes(DemeSettings *demeParameters,
						      thrust::host_vector<int> &everybodys_deme,
						      thrust::host_vector<int> &kids_per_mom,
						      thrust::host_vector<int> &current_deme_sizes,
						      thrust::host_vector<int> &maximum_deme_sizes)
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


void EggsNeonates::inherit_genotypes(thrust::host_vector<float> &probability_individuals_become_mothers,
				 thrust::host_vector<float> &probability_individuals_become_fathers)
	{
	get_maternally_derived_genotype(probability_individuals_become_mothers, species->mgenotype, species->fgenotype);
	get_paternally_derived_genotype(probability_individuals_become_fathers, species->mgenotype, species->fgenotype);
	mutate(species->fgenotype, species->mgenotype);
	}


void EggsNeonates::get_maternally_derived_genotype(thrust::host_vector<float> &probability_individuals_become_mothers,
					     thrust::host_vector<float> *&mgenotype,
					     thrust::host_vector<float> *&fgenotype)
	{
	/*
	* Calculate the offspring's genotypes at the maternally inherited loci. This function will determine who the mother is, and generate a haploid gamete from the mother that will be grafted on to the offspring.
	*/
	mating_ThrustProbTable_demes at;
	thrust::host_vector<int> mother_index(Total_Number_of_Neonates);
	thrust::host_vector<float> rand(Total_Number_of_Neonates);
	thrust::host_vector<int> parity(Total_Number_of_Neonates);

/*
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine mothers.
*/
	at.setup(probability_individuals_become_mothers.begin(), probability_individuals_become_mothers.begin() + previous_pop_size);
	for (int i=0; i < rand.size(); i++)
		{
		rand[i] = gsl_rng_uniform(gen);
		}	
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
		for (int j=0; j < rand.size(); j++)
			{
			rand[j] = gsl_rng_uniform(gen);
			}
		recombine(rand, mother_index, parity, fgenotype, mgenotype, fgenotype[i], i);
		}
	}

void EggsNeonates::get_paternally_derived_genotype(thrust::host_vector<float> &probability_individuals_become_fathers,
					     thrust::host_vector<float> *&mgenotype,
					     thrust::host_vector<float> *&fgenotype)
	{
	/*
	* Calculate the offspring's genotypes at the paternally inherited loci. This function will determine who the mother is, and generate a haploid gamete from the mother that will be grafted on to the offspring.
	*/
	mating_ThrustProbTable_demes at;
	thrust::host_vector<int> father_index(Total_Number_of_Neonates);
	thrust::host_vector<float> rand(Total_Number_of_Neonates);
	thrust::host_vector<int> parity(Total_Number_of_Neonates);
	
	/*
	Directly use phenotype[1] as the reproductive probability.
	Reproductive probability = phenotype[1].
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine fathers.
	*/
	at.setup(probability_individuals_become_fathers.begin(), probability_individuals_become_fathers.begin() + previous_pop_size);

	for (int i=0; i < rand.size(); i++)
		{
		rand[i] = gsl_rng_uniform(gen);
		}

	at.determine_key_offsets( Num_Demes, species->deme_sizes );	
	at.adjust_randoms(rand.begin(), rand.end(), kids_deme.begin(), kids_deme.end());

	at.draw(rand.begin(), rand.end(), father_index.begin());
	thrust::copy(father_index.begin(), father_index.end(), species->paternal_id.begin() + previous_pop_size);


	//Reset parity to zeroes
	thrust::fill(parity.begin(), parity.end(), 0);

	//Recombination for mgenotype
	for (int i = 0 ; i < nloci ; i++) 
		{
		for (int j=0; j < rand.size(); j++)
			{
			rand[j] = gsl_rng_uniform(gen);
			}
		recombine(rand, father_index, parity, fgenotype, mgenotype, mgenotype[i], i);
		}
	}

void EggsNeonates::mutate(thrust::host_vector<float> *&mgenotype, thrust::host_vector<float> *&fgenotype)
	{
	/*
	* Determine, for each offspring's locus, whether there will be a mutation at that locus, and if so what the magnitude of that mutation will be and how that changes the offspring's allelic value. The current behavior assumes mutations are gaussian about the parental allelic value, and that the mutation parameters (mutation rate and sd of gaussian) vary by deme. Possible expansions include allowing the mutation rate itself to be an individual-specific phenotype, or alternative mutational models (e.g., point mutations that have a categorical rather than quantitative effect.
	*/
	thrust::host_vector<float> mutation_size(2*Total_Number_of_Neonates); 

	thrust::copy(mutation_magnitude.begin(), mutation_magnitude.begin() + Total_Number_of_Neonates, mutation_size.begin());
	thrust::copy(mutation_magnitude.begin(), mutation_magnitude.begin() + Total_Number_of_Neonates, mutation_size.begin() + Total_Number_of_Neonates);

	thrust::host_vector<float> mutation(2*Total_Number_of_Neonates);

	// no mutation at the original sex determining locus
	for (int i = 1 ; i < nloci ; i++) 
		{
		thrust::gather(kids_deme.begin(), kids_deme.begin() + Total_Number_of_Neonates, species->demeParameters->GeneticArchitecture->get_mutation_magnitudes_ptr(i), mutation_magnitude.begin());
		thrust::gather(kids_deme.begin(), kids_deme.begin() + Total_Number_of_Neonates, species->demeParameters->GeneticArchitecture->get_mutation_rates_ptr(i), mutation_rate.begin());
		
		for (int j=0; j < Total_Number_of_Neonates; j++)
			{
			int mutation_event = gsl_ran_bernoulli(gen, mutation_rate[j]); 
			if (mutation_event == 1)
				{
				mgenotype[i][j] = mgenotype[i][j] + gsl_ran_gaussian(gen, mutation_magnitude[j]);
				}
		
			mutation_event = gsl_ran_bernoulli(gen, mutation_magnitude[j]);

			if (mutation_event == 1)
				{
				fgenotype[i][j] = fgenotype[i][j] + gsl_ran_gaussian(gen, mutation_magnitude[j]);
				}
			}
		}
	}


void EggsNeonates::recombine(thrust::host_vector<float> &rand,
			     thrust::host_vector<int> &parent,
			     thrust::host_vector<int> &parity,
			     thrust::host_vector<float> *&parents_fgenotype,
			     thrust::host_vector<float> *&parents_mgenotype,
			     thrust::host_vector<float> &kids_genotype,
			     int locus_ID)
	{
/*
* A method that implements recombination. Note that in many cases, the arguments for either (parents_fgenotype or parents_mgenotype) and kids_genotype will be the same vector. This function increments the kids_genotype to begin at the (previous_pop_size)th index. This should probably be refactored so that there will be a kids_mgenotype and kids_fgenotype vector that is local to the neonates class, and that later gets copied into the inds class, i.e. something like: 
\code{.cpp}
	thrust::copy(kids_fgenotype[i].begin(), kids_fgenotype[i].end(), fgenotype[i].begin() + size);
\endcode
*/
	//Set up recombination functor.
	float *fgenotype_ptr = &parents_fgenotype[locus_ID][0];
	float *mgenotype_ptr = &parents_mgenotype[locus_ID][0];
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
