#include "Fish_genotype_phenotype_maps.h"

GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
	{
	if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE"]) 
		{
		return new fecundity_phenotype(species, phenotype_index, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE"])
		{
		return new mortality_phenotype(species, phenotype_index, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["IRREVERSIBLE_MASS_PHENOTYPE"])
		{
		return new irreversible_mass_at_birth(species, phenotype_index, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["REVERSIBLE_MASS_PHENOTYPE"])
		{
		return new reversible_mass_at_birth(species, phenotype_index, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["RESOURCE_LIMITATION_PHENOTYPE"])
		{
		return new satiation_at_birth(species, phenotype_index, index_case, num_kids);
		}

	if (phenotype_index == species->demeParameters->species_specific_values["EGGSIZE_PHENOTYPE"])
		{
		return new eggsize_at_birth(species, phenotype_index, index_case, num_kids);
		}
	}


