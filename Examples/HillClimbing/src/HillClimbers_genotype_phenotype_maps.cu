#include "HillClimbers_genotype_phenotype_maps.h"
#include <math/random_variables_functions.h>

#include <thrust/sequence.h>

#define MAX_MORTALITY 0.975
#define BASELINE_FECUNDITY 2

void fecundity_genotype_phenotype_map::calculate_phenotype(inds *species)
	{
	int deme_offset = 0;
	float coefficient_0 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF0") + deme_offset);
	float coefficient_1 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF1") + deme_offset);
	float coefficient_2 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF2") + deme_offset);
	float coefficient_3 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF3") + deme_offset);
	float coefficient_4 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF4") + deme_offset);

	for (int i=0; i < num_kids; i++)
		{
		if (deme_offset !=  species->deme[index_case + i])
			{
			deme_offset = species->deme[index_case + i];
			coefficient_0 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF0") + deme_offset);
			coefficient_1 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF1") + deme_offset);
			coefficient_2 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF2") + deme_offset);
			coefficient_3 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF3") + deme_offset);
			coefficient_4 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF4") + deme_offset);
			}
		species->phenotype[phenotype_index][index_case + i] = 0.5*coefficient_0*(species->fgenotype[0][index_case + i] + species->mgenotype[0][index_case + i]) 
								 + 0.5*coefficient_1*(species->fgenotype[1][index_case + i] + species->mgenotype[1][index_case + i]) 
								 + 0.5*coefficient_2*(species->fgenotype[2][index_case + i] + species->mgenotype[2][index_case + i]) 
								 + 0.5*coefficient_3*(species->fgenotype[3][index_case + i] + species->mgenotype[3][index_case + i]) 
								 + 0.5*coefficient_4*(species->fgenotype[4][index_case + i] + species->mgenotype[4][index_case + i]);

		species->phenotype[phenotype_index][index_case + i] = fabs(species->phenotype[phenotype_index][index_case + i]) + BASELINE_FECUNDITY;
		}
	}

void mortality_genotype_phenotype_map::calculate_phenotype(inds *species)
	{	
	int deme_offset = 0;
	float coefficient_0 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF0") + deme_offset);
	float coefficient_1 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF1") + deme_offset);
	float coefficient_2 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF2") + deme_offset);
	float coefficient_3 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF3") + deme_offset);
	float coefficient_4 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF4") + deme_offset);

	
	for (int i=0; i < num_kids; i++)
		{
		if (deme_offset !=  species->deme[index_case + i])
			{
			deme_offset = species->deme[index_case + i];
			coefficient_0 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF0") + deme_offset);
			coefficient_1 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF1") + deme_offset);
			coefficient_2 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF2") + deme_offset);
			coefficient_3 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF3") + deme_offset);
			coefficient_4 = *(Parameters->get_vector_ptr("GENPHEN_MAP_COEF4") + deme_offset);
			}
		species->phenotype[phenotype_index][index_case + i] = 0.5*coefficient_0*(species->fgenotype[5][index_case + i] + species->mgenotype[5][index_case + i]) 
								 + 0.5*coefficient_1*(species->fgenotype[6][index_case + i] + species->mgenotype[6][index_case + i]) 
								 + 0.5*coefficient_2*(species->fgenotype[7][index_case + i] + species->mgenotype[7][index_case + i]) 
								 + 0.5*coefficient_3*(species->fgenotype[8][index_case + i] + species->mgenotype[8][index_case + i]) 
								 + 0.5*coefficient_4*(species->fgenotype[9][index_case + i] + species->mgenotype[9][index_case + i]);

		// Make sure the resulting value is biologically meaningful:
		species->phenotype[phenotype_index][index_case + i] = exp(-1*fabs(species->phenotype[phenotype_index][index_case + i]));

		if (species->phenotype[phenotype_index][index_case + i] > MAX_MORTALITY)
			{
			species->phenotype[phenotype_index][index_case + i] = MAX_MORTALITY;
			}
		}
	}
 
