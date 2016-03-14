#include "myResource_genotype_phenotype_maps.h"
#include "myConsumer_genotype_phenotype_maps.h"

GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
	{
	if (species->species_ID == 0)
		{
		if (phenotype_index == species->demeParameters->species_specific_values["CONSUMER_ATTACK_PHENOTYPE_INDEX"])
			{
			return new consumer_attack_phenotype( species, index_case, num_kids );
			}
		if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"])
			{
			return new consumer_mortality_phenotype( species, index_case, num_kids );
			}
		if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"])
			{
			return new consumer_fecundity_phenotype( species, index_case, num_kids );
			}
		}
	if (species->species_ID == 1)
		{
		if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"])
			{
			return new resource_fecundity_phenotype( species, index_case, num_kids );
			}
		if (phenotype_index == species->demeParameters->species_specific_values["RESOURCE_DEFENSE_PHENOTYPE_INDEX"])
			{
			return new resource_defense_phenotype( species, index_case, num_kids );
			}
		if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"])
			{
			return new resource_mortality_phenotype( species, index_case, num_kids );
			}
		if (phenotype_index == species->demeParameters->species_specific_values["RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX"])
			{
			return new resource_competition_phenotype( species, index_case, num_kids );
			}
		}
	}


