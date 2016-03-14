#include "HillClimbers_genotype_phenotype_maps.h"

GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
	{
	if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]) 
		{
		return new fecundity_genotype_phenotype_map(species, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"])
		{
		return new mortality_genotype_phenotype_map(species, index_case, num_kids);
		}
	}


