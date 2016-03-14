#include <species/add_kids/assortative_mating_neonates_class.h>

/* One problem with a separate neonates class is that all this stuff has to get copied back into example species. It might just not be worth it depending on the performance cost. */

#define NLOCI 50

float recomb_array_assort[NLOCI] = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,0.5f, 0.5f, 0.5f, 0.5f, 0.5f};


Assortative_mating_neonates::Assortative_mating_neonates(thrust::device_vector<int> &pair_populations, DemeSettings *subpopParameters, thrust::device_vector<int> &everybodys_deme, thrust::device_vector<int> &kids_per_mom,  thrust::device_vector<int> &current_deme_sizes, thrust::device_vector<int> &maximum_deme_sizes, int N_alive_inds,  int num_loci, int nPhen) : EggsNeonates(subpopParameters, everybodys_deme, kids_per_mom, current_deme_sizes, maximum_deme_sizes, N_alive_inds, num_loci, nPhen) 
	{

	pairs_per_deme.resize(Num_Subpopulations);
	
	// Determine how many pairs are in each subpopulation
	thrust::counting_iterator<int> search_begin(0);
	thrust::device_vector<int> temp_Subpop_sizes;

	temp_Subpop_sizes.resize(Num_Subpopulations);

	thrust::upper_bound(pair_populations.begin(), pair_populations.end(),
                      search_begin, search_begin + Num_Subpopulations,
                      temp_Subpop_sizes.begin());

	thrust::adjacent_difference(temp_Subpop_sizes.begin(), temp_Subpop_sizes.end(),
                              pairs_per_deme.begin());
		
	}

void Assortative_mating_neonates::inherit_genotypes_by_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
				thrust::device_vector<int> &fathers_list,
				thrust::device_vector<int> &mothers_list,
				thrust::device_vector<float> *&fgenotype,
				thrust::device_vector<float> *&mgenotype,
				curandGenerator_t generator)
	{

	mothers_chosen.resize(Total_Number_of_Neonates);
	fathers_chosen.resize(Total_Number_of_Neonates);

	get_mating_pair(probability_pair_becomes_parents, fathers_list, mothers_list, generator);

	get_maternally_derived_genotype_deterministic(mothers_chosen, mgenotype, fgenotype, generator);

	get_paternally_derived_genotype_deterministic(fathers_chosen, mgenotype, fgenotype, generator);

	mutate(generator, fgenotype, mgenotype);
	}

void Assortative_mating_neonates::get_mating_pair(thrust::device_vector<float> &probability_pair_becomes_parents,
						   thrust::device_vector<int> &fathers_list,
						   thrust::device_vector<int> &mothers_list,
						   curandGenerator_t generator
						  )
	{


	mating_subpopThrustProbTable at;
	thrust::device_vector<int> pair_index(Total_Number_of_Neonates);
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
/*
	Feed reproductive probablity into the setup of the alias table.
	Draw from the alias table to determine mothers.
*/
	at.setup(probability_pair_becomes_parents.begin(), probability_pair_becomes_parents.end());

	curandGenerateUniform(generator, rand_ptr, Total_Number_of_Neonates);
 
	at.determine_key_offsets( Num_Subpopulations, pairs_per_deme );
 
	at.adjust_randoms(rand.begin(), rand.end(), kids_pop.begin(), kids_pop.end());
 
	at.draw(rand.begin(), rand.end(), pair_index.begin());

 		
	thrust::gather(pair_index.begin(), pair_index.end(), fathers_list.begin(), fathers_chosen.begin());

 
	thrust::gather(pair_index.begin(), pair_index.end(), mothers_list.begin(), mothers_chosen.begin());

	}

/* use the functions blah_blah_deterministic() if the parents have already been chosen and you just need to copy genotypes */
	
void Assortative_mating_neonates::get_maternally_derived_genotype_deterministic(thrust::device_vector<int> &mother_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype,
					     curandGenerator_t generator)
	{
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(generator, rand_ptr, Total_Number_of_Neonates);
	
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	thrust::fill(parity.begin(), parity.end(), 0);

	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(generator, rand_ptr, Total_Number_of_Neonates);
                recombine(rand, mother_index, parity, fgenotype, mgenotype, fgenotype[i], i);
		}
	}

void Assortative_mating_neonates::get_paternally_derived_genotype_deterministic(thrust::device_vector<int> &father_index,
					     thrust::device_vector<float> *&mgenotype,
					     thrust::device_vector<float> *&fgenotype,
					     curandGenerator_t generator)
	{
	thrust::device_vector<float> rand(Total_Number_of_Neonates);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(generator, rand_ptr, Total_Number_of_Neonates);
	
	//Reset parity to zeroes
	thrust::device_vector<int> parity(Total_Number_of_Neonates);
	thrust::fill(parity.begin(), parity.end(), 0);

	//Recombination for mgenotype
	for (int i = 0 ; i < nloci ; i++) 
		{
		curandGenerateUniform(generator, rand_ptr,Total_Number_of_Neonates);
		recombine(rand, father_index, parity, fgenotype, mgenotype, mgenotype[i], i);
		}
	}

void Assortative_mating_neonates::record_parents(thrust::device_vector<int> &maternal_id, thrust::device_vector<int> &paternal_id, thrust::device_vector<int> &ids)
	{
	thrust::gather(mothers_chosen.begin(), mothers_chosen.begin() + Total_Number_of_Neonates, ids.begin(), maternal_id.begin() + current_pop_size);
	thrust::gather(fathers_chosen.begin(), fathers_chosen.begin() + Total_Number_of_Neonates, ids.begin(), paternal_id.begin() + current_pop_size);
	}
