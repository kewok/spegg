#include <species/add_kids/genotype_phenotype_map.h>
#include <species/inds_stochastic.h>

class fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		fecundity_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};

class mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		mortality_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};

class irreversible_mass_at_birth : public GenotypePhenotypeMap
	{
	public:
		irreversible_mass_at_birth(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};


class reversible_mass_at_birth : public GenotypePhenotypeMap
	{
	public:
		reversible_mass_at_birth(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};


class satiation_at_birth : public GenotypePhenotypeMap
	{
	public:
		satiation_at_birth(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};


class eggsize_at_birth : public GenotypePhenotypeMap
	{
	public:
		eggsize_at_birth(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			gen = (static_cast<inds_stochastic*> (species))->gen;
			}
	
		void calculate_phenotype(inds *species);
	};


/* define functors, if you have any */

struct offspring_size_calculator
{
float *intercept;

float *fgen1;
float *fgen2;
float *fgen3;
float *fgen4;
float *fgen5;
float *fgen6;
float *fgen7;

float *mgen1;
float *mgen2;
float *mgen3;
float *mgen4;
float *mgen5;
float *mgen6;
float *mgen7;

float *allelic_effect1;
float *allelic_effect2;
float *allelic_effect3;
float *allelic_effect4;
float *allelic_effect5;
float *allelic_effect6;
float *allelic_effect7;

float *maternal_effects;

offspring_size_calculator(float *intercept_ptr, float *fgen_ptr1,float *fgen_ptr2,float *fgen_ptr3,float *fgen_ptr4,float *fgen_ptr5, float *fgen_ptr6, float *fgen_ptr7, float *mgen_ptr1, float *mgen_ptr2, float *mgen_ptr3, float *mgen_ptr4, float *mgen_ptr5, float *mgen_ptr6, float *mgen_ptr7, float *allelic_effect_ptr1, float *allelic_effect_ptr2, float *allelic_effect_ptr3, float *allelic_effect_ptr4, float *allelic_effect_ptr5, float *allelic_effect_ptr6, float *allelic_effect_ptr7, float *maternal_effects_ptr) : intercept(intercept_ptr), fgen1(fgen_ptr1), fgen2(fgen_ptr2), fgen3(fgen_ptr3), fgen4(fgen_ptr4), fgen5(fgen_ptr5), fgen6(fgen_ptr6), fgen7(fgen_ptr7), mgen1(mgen_ptr1), mgen2(mgen_ptr2), mgen3(mgen_ptr3), mgen4(mgen_ptr4), mgen5(mgen_ptr5), mgen6(mgen_ptr6), mgen7(mgen_ptr7), allelic_effect1(allelic_effect_ptr1), allelic_effect2(allelic_effect_ptr2), allelic_effect3(allelic_effect_ptr3), allelic_effect4(allelic_effect_ptr4), allelic_effect5(allelic_effect_ptr5), allelic_effect6(allelic_effect_ptr6), allelic_effect7(allelic_effect_ptr7), maternal_effects(maternal_effects_ptr)
	{};

	/* 
		Elements in the tuple.

		----------------------
		0: individual index
		1: random noise
		2: phenotype
		3: individuals deme
	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		int ind_deme = thrust::get<3>(t);
		thrust::get<2>(t) = intercept[ind_deme] + 
				allelic_effect1[ind_deme]*0.5*(fgen1[thrust::get<0>(t)] + mgen1[thrust::get<0>(t)]) +
				allelic_effect2[ind_deme]*0.5*(fgen2[thrust::get<0>(t)] + mgen2[thrust::get<0>(t)]) + 
				allelic_effect3[ind_deme]*0.5*(fgen3[thrust::get<0>(t)] + mgen3[thrust::get<0>(t)]) + 
				allelic_effect4[ind_deme]*0.5*(fgen4[thrust::get<0>(t)] + mgen4[thrust::get<0>(t)]) + 
				allelic_effect5[ind_deme]*0.5*(fgen5[thrust::get<0>(t)] + mgen5[thrust::get<0>(t)]) + 
				allelic_effect6[ind_deme]*0.5*(fgen6[thrust::get<0>(t)] + mgen6[thrust::get<0>(t)]) + 
				allelic_effect7[ind_deme]*0.5*(fgen7[thrust::get<0>(t)] + mgen7[thrust::get<0>(t)]) + 		thrust::get<1>(t);

		/* transform to be no smaller than 0.001 g */
		if (thrust::get<2>(t) < 0.001)
			{
			thrust::get<2>(t) = 0.001;
			}
		}
	
};

