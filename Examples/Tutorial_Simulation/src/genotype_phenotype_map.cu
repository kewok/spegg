#include "Penguins_genotype_phenotype_maps.h"

GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
    {
    if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]) 
	{       
	return new fecundity_phenotype(species, phenotype_index, index_case, num_kids);
	}

    if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"])
	{        
	return new mortality_phenotype(species, phenotype_index, index_case, num_kids);
	}

   if (phenotype_index == species->demeParameters->species_specific_values["CROWN_COLOR_INDEX"])
	{
        return new crown_color_phenotype(species, phenotype_index, index_case, num_kids);
	}
    }
