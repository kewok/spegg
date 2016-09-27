#ifndef GENOTYPE_PHENOTYPE_MAP_PARAMETERS_H
#define GENOTYPE_PHENOTYPE_MAP_PARAMETERS_H

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <algorithm>

using namespace libconfig;

class GenotypePhenotypeMapParameters
	{
	friend class DemeGeneticsSettings;

	GenotypePhenotypeMapParameters(const char *filename, int species_ID, int phenotype_index, std::vector<std::string> &loci_names);

	public: 
		thrust::device_vector<float> *deme_specific_parameters;

		thrust::device_ptr<float> get_vector_ptr(const char *parameter_name);
	
	protected:
		int Number_of_Parameters; 
		int phenotype_index;

		std::vector<std::string> Names_of_Genotype_Phenotype_Map_Parameters;
		int Number_of_Demes;

		void read_in_data(const char *filename, int species_ID);

		void specify_parameter_index();

		std::map<std::string, int> parameter_index;
	};

#endif
