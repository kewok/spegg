#ifndef GENETIC_DEME_SPECIFIC_DATA_CLASS_H
#define GENETIC_DEME_SPECIFIC_DATA_CLASS_H

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <map>
#include <thrust/functional.h>

#include <species/add_kids/genotype_phenotype_map_parameters.h>

using namespace libconfig;

class DemeGeneticsSettings
{
public: 
	DemeGeneticsSettings(const char *filename, int species_ID);

	thrust::host_vector<float> recombination_rates;
	GenotypePhenotypeMapParameters **phen_gen_map_parm;

	float* get_vector_ptr(int locus_number);

	std::vector<std::string> phenotype_names;

	float *get_mutation_rates_ptr(int locus_index);
	float *get_mutation_magnitudes_ptr(int locus_index);

	int Number_of_Loci;
	int Number_of_Phenotypes;

protected:
	void read_in_data(const char *filename, int species_ID);
	int Number_of_Demes;

	thrust::host_vector<float> *deme_specific_mutation_rates;
	thrust::host_vector<float> *deme_specific_mutation_magnitudes;

	std::vector<std::string> loci_names;
};

#endif
