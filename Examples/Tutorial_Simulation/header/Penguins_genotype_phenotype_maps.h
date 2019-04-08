#include <species/add_kids/genotype_phenotype_map.h>

class fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		fecundity_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};
	
		void calculate_phenotype(inds *species);
	};

struct fecundity_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;
	float *map_coefficient_1;

	fecundity_calculator(float* map_cons, float* map_coef0, float* map_coef1) : map_constant(map_cons), map_coefficient_0(map_coef0), map_coefficient_1(map_coef1)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's maternally inherited allelic value at their locus 1
	4: individual's paternally inherited allelic value at their locus 1
	5: individual's fecundity phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<5>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 + map_coefficient_1[ind_deme] * (thrust::get<3>(t) + thrust::get<4>(t) )/2;
		}   
	};

class mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		mortality_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};

		void calculate_phenotype(inds *species);
	};


struct mortality_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;

	mortality_calculator(float* map_cons, float* map_coef0) : map_constant(map_cons), map_coefficient_0(map_coef0)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's mortality phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<3>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 ;
		}   
	};


class crown_color_phenotype : public GenotypePhenotypeMap
	{
	public:
		int CROWN_COLOR_INDEX;

		crown_color_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			this->CROWN_COLOR_INDEX = this->phenotype_index;
			};

		void calculate_phenotype(inds *species);
	};


struct crown_color_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;

	crown_color_calculator(float* map_cons, float* map_coef0) : map_constant(map_cons), map_coefficient_0(map_coef0)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's crown color phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<3>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 ;
		}   
	};

