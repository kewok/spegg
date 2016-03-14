#include <thrust/gather.h>
#include <thrust/sequence.h>

#include "myResource_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

void resource_fecundity_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[FECUNDITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void resource_mortality_phenotype::calculate_phenotype(inds *species)
	{
	thrust::gather(species->deme.begin() + index_case, species->deme.begin() + index_case + num_kids, Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT"), species->phenotype[MORTALITY_PHENOTYPE_INDEX].begin() + index_case);
	}

void resource_defense_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = (&species->fgenotype[0][0]);
	float *fgen2 = (&species->fgenotype[1][0]);
	float *fgen3 = (&species->fgenotype[2][0]);
	float *fgen9 = (&species->fgenotype[8][0]);
	float *fgen10 = (&species->fgenotype[9][0]);

	float *mgen1 = (&species->mgenotype[0][0]);
	float *mgen2 = (&species->mgenotype[1][0]);
	float *mgen3 = (&species->mgenotype[2][0]);
	float *mgen9 = (&species->mgenotype[8][0]);
	float *mgen10 = (&species->mgenotype[9][0]);

// wrap the allelic effects:
	float *locus_effect1 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *locus_effect2 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *locus_effect3 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *locus_effect9 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF9")[0]);
	float *locus_effect10 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF10")[0]);

	int deme_offset = 0;
	float COEF_1 = *(locus_effect1 + deme_offset);
	float COEF_2 = *(locus_effect2 + deme_offset);
	float COEF_3 = *(locus_effect3 + deme_offset);
	float COEF_4 = *(locus_effect9 + deme_offset);
	float COEF_5 = *(locus_effect10 + deme_offset);

	for (int i=index_case; i < index_case + num_kids; i++)
		{
		if (species->deme[i] != deme_offset)
			{
			deme_offset = species->deme[i];

			COEF_1 = *(locus_effect1 + deme_offset);
			COEF_2 = *(locus_effect2 + deme_offset);
			COEF_3 = *(locus_effect3 + deme_offset);
			COEF_4 = *(locus_effect9 + deme_offset);
			COEF_5 = *(locus_effect10 + deme_offset);
			}
		species->phenotype[RESOURCE_DEFENSE_PHENOTYPE_INDEX][i] = COEF_1*(species->fgenotype[0][i] + species->mgenotype[0][i]) + COEF_2*(species->fgenotype[1][i] + species->mgenotype[1][i]) + COEF_3*(species->fgenotype[2][i] + species->mgenotype[2][i]) + COEF_4*(species->fgenotype[3][i] + species->mgenotype[3][i]) + COEF_5*(species->fgenotype[4][i] + species->mgenotype[4][i]);	
		species->phenotype[RESOURCE_DEFENSE_PHENOTYPE_INDEX][i] = species->phenotype[RESOURCE_DEFENSE_PHENOTYPE_INDEX][i]*0.5;
		}
	}

void resource_competition_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the genotypes:
	float *fgen1 = (&species->fgenotype[0][0]);
	float *fgen2 = (&species->fgenotype[1][0]);
	float *fgen3 = (&species->fgenotype[2][0]);
	float *fgen8 = (&species->fgenotype[7][0]);
	float *fgen10 = (&species->fgenotype[9][0]);

	float *mgen1 = (&species->mgenotype[0][0]);
	float *mgen2 = (&species->mgenotype[1][0]);
	float *mgen3 = (&species->mgenotype[2][0]);
	float *mgen8 = (&species->mgenotype[7][0]);
	float *mgen10 = (&species->mgenotype[9][0]);

// wrap the allelic effects:
	float *locus_effect1 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);
	float *locus_effect2 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF2")[0]);
	float *locus_effect3 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF3")[0]);
	float *locus_effect8 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF8")[0]);
	float *locus_effect10 = (&Parameters->get_vector_ptr("GENPHEN_MAP_COEF10")[0]);

	int deme_offset = 0;
	float COEF_1 = *(locus_effect1 + deme_offset);
	float COEF_2 = *(locus_effect2 + deme_offset);
	float COEF_3 = *(locus_effect3 + deme_offset);
	float COEF_4 = *(locus_effect8 + deme_offset);
	float COEF_5 = *(locus_effect10 + deme_offset);

	for (int i=index_case; i < index_case + num_kids; i++)
		{
		if (species->deme[i] != deme_offset)
			{
			deme_offset = species->deme[i];

			COEF_1 = *(locus_effect1 + deme_offset);
			COEF_2 = *(locus_effect2 + deme_offset);
			COEF_3 = *(locus_effect3 + deme_offset);
			COEF_4 = *(locus_effect8 + deme_offset);
			COEF_5 = *(locus_effect10 + deme_offset);
			}
		species->phenotype[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX][i] = COEF_1*(species->fgenotype[0][i] + species->mgenotype[0][i]) + COEF_2*(species->fgenotype[1][i] + species->mgenotype[1][i]) + COEF_3*(species->fgenotype[2][i] + species->mgenotype[2][i]) + COEF_4*(species->fgenotype[3][i] + species->mgenotype[3][i]) + COEF_5*(species->fgenotype[4][i] + species->mgenotype[4][i]);	
		species->phenotype[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX][i] = species->phenotype[RESOURCE_COMPETITIVE_ABILITY_PHENOTYPE_INDEX][i]*0.5;
		}
	}
