#include <species/add_kids/genotype_phenotype_map.h>

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
