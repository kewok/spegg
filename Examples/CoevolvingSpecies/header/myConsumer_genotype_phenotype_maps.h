#ifndef MY_CONSUMER_GENOTYPE_PHENOTYPE_MAP_H
#define MY_CONSUMER_GENOTYPE_PHENOTYPE_MAP_H

#include <species/add_kids/genotype_phenotype_map.h>

class consumer_fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		consumer_fecundity_phenotype(inds *species, int index_case, int num_kids)
			{
			this->FECUNDITY_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[FECUNDITY_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int FECUNDITY_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};

class consumer_mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		consumer_mortality_phenotype(inds *species, int index_case, int num_kids)
			{
			this->MORTALITY_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[MORTALITY_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int MORTALITY_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};


class consumer_attack_phenotype : public GenotypePhenotypeMap
	{
	public:
		consumer_attack_phenotype(inds *species, int index_case, int num_kids)
			{
			this->CONSUMER_ATTACK_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["CONSUMER_ATTACK_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[CONSUMER_ATTACK_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int CONSUMER_ATTACK_PHENOTYPE_INDEX;
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
