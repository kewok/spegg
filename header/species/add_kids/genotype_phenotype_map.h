#ifndef GENOTYPE_PHENOTYPE_MAP
#define GENOTYPE_PHENOTYPE_MAP

#include <species/inds.h>
#include <species/add_kids/genotype_phenotype_map_parameters.h>

//virtual interface

class GenotypePhenotypeMap
	{
	public:
/* use a factory method */
		static GenotypePhenotypeMap *create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids);

		GenotypePhenotypeMap(inds *species, int phenotype_index, int index_case, int num_kids)
			{
			this->phenotype_index = phenotype_index;
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[phenotype_index];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
	
		virtual void calculate_phenotype(inds *species)=0;

	protected:
		GenotypePhenotypeMapParameters *Parameters;

		int phenotype_index;
		int index_case;
		int num_kids;

		curandGenerator_t gen;		
	};

#endif
