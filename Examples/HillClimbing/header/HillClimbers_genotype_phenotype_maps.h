#include <species/add_kids/genotype_phenotype_map.h>

#define BASELINE_FECUNDITY 2

class fecundity_genotype_phenotype_map : public GenotypePhenotypeMap
	{
	public:
		fecundity_genotype_phenotype_map(inds *species, int index_case, int num_kids)
			{
			this->index_case = index_case;
			this->num_kids = num_kids;

			this->phenotype_index = species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[phenotype_index];
			}

		void calculate_phenotype(inds *species);
	};


struct fecundity_calculator
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

fecundity_calculator(float *fgen_ptr1,float *fgen_ptr2,float *fgen_ptr3,float *fgen_ptr4,float *fgen_ptr5, float *mgen_ptr1, float *mgen_ptr2, float *mgen_ptr3, float *mgen_ptr4, float *mgen_ptr5, float *allelic_effect_ptr1, float *allelic_effect_ptr2, float *allelic_effect_ptr3, float *allelic_effect_ptr4, float *allelic_effect_ptr5) : fgen1(fgen_ptr1), fgen2(fgen_ptr2), fgen3(fgen_ptr3), fgen4(fgen_ptr4), fgen5(fgen_ptr5), mgen1(mgen_ptr1), mgen2(mgen_ptr2), mgen3(mgen_ptr3), mgen4(mgen_ptr4), mgen5(mgen_ptr5), allelic_effect1(allelic_effect_ptr1), allelic_effect2(allelic_effect_ptr2), allelic_effect3(allelic_effect_ptr3), 
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

		/* transform to give numbers >= 1*/
		thrust::get<1>(t) = fabs(thrust::get<1>(t)) + BASELINE_FECUNDITY;
		}
};


class mortality_genotype_phenotype_map : public GenotypePhenotypeMap
	{
	public:
		mortality_genotype_phenotype_map(inds *species, int index_case, int num_kids )
			{
			this->phenotype_index = species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"];
			this->Parameters = species->demeParameters->GeneticArchitecture -> phen_gen_map_parm[phenotype_index];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}

		void calculate_phenotype(inds *species);
	};


struct mortality_calculator
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

mortality_calculator(float *fgen_ptr1,float *fgen_ptr2,float *fgen_ptr3,float *fgen_ptr4,float *fgen_ptr5, float *mgen_ptr1, float *mgen_ptr2, float *mgen_ptr3, float *mgen_ptr4, float *mgen_ptr5, float *allelic_effect_ptr1, float *allelic_effect_ptr2, float *allelic_effect_ptr3, float *allelic_effect_ptr4, float *allelic_effect_ptr5) : fgen1(fgen_ptr1), fgen2(fgen_ptr2), fgen3(fgen_ptr3), fgen4(fgen_ptr4), fgen5(fgen_ptr5), mgen1(mgen_ptr1), mgen2(mgen_ptr2), mgen3(mgen_ptr3), mgen4(mgen_ptr4), mgen5(mgen_ptr5), allelic_effect1(allelic_effect_ptr1), allelic_effect2(allelic_effect_ptr2), allelic_effect3(allelic_effect_ptr3), 
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

		// Make sure the resulting value is biologically meaningful.
		thrust::get<1>(t) = exp(-1*fabs(thrust::get<1>(t)));
		}
};
