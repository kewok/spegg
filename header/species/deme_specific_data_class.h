#ifndef DEME_SPECIFIC_DATA_CLASS_H
#define DEME_SPECIFIC_DATA_CLASS_H

#include <species/add_kids/genetic_deme_specific_data.h>

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <libconfig.h++>
#include <string>
#include <curand.h>
#include <map>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <algorithm>

using namespace libconfig;

class DemeSettings
	{
	public: 
		DemeSettings(const char *filename, int species_ID);

		thrust::device_vector<float> *deme_wide_parameters;

		std::map<std::string, int> parameter_index;
		std::map<std::string, float> species_specific_values;

		DemeGeneticsSettings *GeneticArchitecture;

		thrust::device_vector<float> *deme_specific_recombination_rates;

		thrust::device_ptr<float> get_vector_ptr(const char *parameter_name);
		thrust::device_ptr<float> get_allelic_effects_ptr(int locus_number);

		bool does_parameter_exist(const char *parameter_name);

		int check_number_of_demes();
	protected:
		void specify_parameter_index();
		void read_in_parameters(const char *filename, int species_ID);
		int Number_of_Demes;
		int Number_of_Parameters;
		int Number_of_Species_Specific_Values;

		std::vector<std::string> parameter_names;	
		std::vector<std::string> species_specific_values_names;
	};

#endif
