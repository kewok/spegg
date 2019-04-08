#ifndef MY_CONSUMER_GENOTYPE_PHENOTYPE_MAP_H
#define MY_CONSUMER_GENOTYPE_PHENOTYPE_MAP_H

#include <species/add_kids/genotype_phenotype_map.h>

class consumer_fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		consumer_fecundity_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};
	
		void calculate_phenotype(inds *species);
	};

class consumer_mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		consumer_mortality_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};
	
		void calculate_phenotype(inds *species);
	};


class consumer_attack_phenotype : public GenotypePhenotypeMap
	{
	public:
		int CONSUMER_ATTACK_PHENOTYPE_INDEX;

		consumer_attack_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			CONSUMER_ATTACK_PHENOTYPE_INDEX = this->phenotype_index;
			};

		void calculate_phenotype(inds *species);
	};

struct consumer_attack_phenotype_calculator
	{
	float *map_coefficient_0;
	float *map_coefficient_1;
	float *map_coefficient_2;

	consumer_attack_phenotype_calculator(float *map_coef0, float *map_coef1, float *map_coef2) : map_coefficient_0(map_coef0), map_coefficient_1(map_coef1), map_coefficient_2(map_coef2)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's maternally inherited allelic value at their locus 1
	4: individual's paternally inherited allelic value at their locus 1
	5: individual's maternally inherited allelic value at their locus 2
	6: individual's paternally inherited allelic value at their locus 2
	7: individual consumer's attack rate phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<7>(t) = map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 + map_coefficient_1[ind_deme] * (thrust::get<3>(t) + thrust::get<4>(t) )/2 + map_coefficient_2[ind_deme] * (thrust::get<5>(t) + thrust::get<6>(t) )/2;
		}   
	};

#endif
