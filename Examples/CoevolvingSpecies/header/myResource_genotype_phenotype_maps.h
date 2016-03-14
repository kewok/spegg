#ifndef MY_RESOURCE_GENOTYPE_PHENOTYPE_MAP_H
#define MY_RESOURCE_GENOTYPE_PHENOTYPE_MAP_H

#include <species/add_kids/genotype_phenotype_map.h>

class resource_mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		resource_mortality_phenotype(inds *species, int index_case, int num_kids)
			{
			this->MORTALITY_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[MORTALITY_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int MORTALITY_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};

class resource_fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		resource_fecundity_phenotype(inds *species, int index_case, int num_kids)
			{
			this->FECUNDITY_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[FECUNDITY_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int FECUNDITY_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};

class resource_defense_phenotype : public GenotypePhenotypeMap
	{
	public:
		resource_defense_phenotype(inds *species, int index_case, int num_kids)
			{
			this->RESOURCE_DEFENSE_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["RESOURCE_DEFENSE_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[RESOURCE_DEFENSE_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int RESOURCE_DEFENSE_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};

class resource_competition_phenotype : public GenotypePhenotypeMap
	{
	public:
		resource_competition_phenotype(inds *species, int index_case, int num_kids)
			{
			this->RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX = species->demeParameters->species_specific_values["RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		int RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX;
		void calculate_phenotype(inds *species);
	};

struct resource_defense_calculator
{
float *fgen1;
float *fgen2;
float *fgen3;
float *fgen4;
float *fgen5;
float *mgen1;
float *mgen2;
float *mgen3;
float *mgen4;
float *mgen5;
float *allelic_effect1;
float *allelic_effect2;
float *allelic_effect3;
float *allelic_effect4;
float *allelic_effect5;

resource_defense_calculator(float *fgen_ptr1,float *fgen_ptr2,float *fgen_ptr3,float *fgen_ptr4,float *fgen_ptr5, float *mgen_ptr1, float *mgen_ptr2, float *mgen_ptr3, float *mgen_ptr4, float *mgen_ptr5, float *allelic_effect_ptr1, float *allelic_effect_ptr2, float *allelic_effect_ptr3, float *allelic_effect_ptr4, float *allelic_effect_ptr5) : fgen1(fgen_ptr1), fgen2(fgen_ptr2), fgen3(fgen_ptr3), fgen4(fgen_ptr4), fgen5(fgen_ptr5), mgen1(mgen_ptr1), mgen2(mgen_ptr2), mgen3(mgen_ptr3), mgen4(mgen_ptr4), mgen5(mgen_ptr5), allelic_effect1(allelic_effect_ptr1), allelic_effect2(allelic_effect_ptr2), allelic_effect3(allelic_effect_ptr3), 
allelic_effect4(allelic_effect_ptr4), allelic_effect5(allelic_effect_ptr5)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: individual index
		1: phenotype
		2: individuals deme
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ind_deme = thrust::get<2>(t);
		thrust::get<1>(t) = allelic_effect1[ind_deme]*0.5*(fgen1[thrust::get<0>(t)] + mgen1[thrust::get<0>(t)]) +
				allelic_effect2[ind_deme]*0.5*(fgen2[thrust::get<0>(t)] + mgen2[thrust::get<0>(t)]) + 
				allelic_effect3[ind_deme]*0.5*(fgen3[thrust::get<0>(t)] + mgen3[thrust::get<0>(t)]) + 
				allelic_effect4[ind_deme]*0.5*(fgen4[thrust::get<0>(t)] + mgen4[thrust::get<0>(t)]) + 
				allelic_effect5[ind_deme]*0.5*(fgen5[thrust::get<0>(t)] + mgen5[thrust::get<0>(t)]);		
		}
};


struct resource_competitive_ability_calculator
{
float *fgen1;
float *fgen2;
float *fgen3;
float *fgen4;
float *fgen5;
float *mgen1;
float *mgen2;
float *mgen3;
float *mgen4;
float *mgen5;
float *allelic_effect1;
float *allelic_effect2;
float *allelic_effect3;
float *allelic_effect4;
float *allelic_effect5;

resource_competitive_ability_calculator(float *fgen_ptr1,float *fgen_ptr2,float *fgen_ptr3,float *fgen_ptr4,float *fgen_ptr5, float *mgen_ptr1, float *mgen_ptr2, float *mgen_ptr3, float *mgen_ptr4, float *mgen_ptr5, float *allelic_effect_ptr1, float *allelic_effect_ptr2, float *allelic_effect_ptr3, float *allelic_effect_ptr4, float *allelic_effect_ptr5) : fgen1(fgen_ptr1), fgen2(fgen_ptr2), fgen3(fgen_ptr3), fgen4(fgen_ptr4), fgen5(fgen_ptr5), mgen1(mgen_ptr1), mgen2(mgen_ptr2), mgen3(mgen_ptr3), mgen4(mgen_ptr4), mgen5(mgen_ptr5), allelic_effect1(allelic_effect_ptr1), allelic_effect2(allelic_effect_ptr2), allelic_effect3(allelic_effect_ptr3), allelic_effect4(allelic_effect_ptr4), allelic_effect5(allelic_effect_ptr5)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: individual index
		1: phenotype
		2: individuals deme
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ind_deme = thrust::get<2>(t);
		thrust::get<1>(t) = allelic_effect1[ind_deme]*0.5*(fgen1[thrust::get<0>(t)] + mgen1[thrust::get<0>(t)]) +
				allelic_effect2[ind_deme]*0.5*(fgen2[thrust::get<0>(t)] + mgen2[thrust::get<0>(t)]) + 
				allelic_effect3[ind_deme]*0.5*(fgen3[thrust::get<0>(t)] + mgen3[thrust::get<0>(t)]) + 
				allelic_effect4[ind_deme]*0.5*(fgen4[thrust::get<0>(t)] + mgen4[thrust::get<0>(t)]) + 
				allelic_effect5[ind_deme]*0.5*(fgen5[thrust::get<0>(t)] + mgen5[thrust::get<0>(t)]);
		}
};

#endif
