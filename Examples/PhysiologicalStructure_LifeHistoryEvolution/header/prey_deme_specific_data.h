#ifndef PREY_DEME_SPECIFIC_DATA_H
#define PREY_DEME_SPECIFIC_DATA_H

#include <libconfig.h++>
#include <string>
#include <map>
#include <thrust/host_vector.h>
#include <thrust/functional.h>

using namespace libconfig;

class PreyDemeSpecificData
	{
	public:
		PreyDemeSpecificData(const char *filename, int prey_index);
		thrust::host_vector<float> *deme_wide_parameters;

		std::map<std::string, int> parameter_index;
		float *get_vector_ptr(const char *parameter_name);

	protected:
		void read_in_prey_parameters(const char *filename, int prey_index);
		void specify_parameter_index();
		int Number_of_Demes;
		int Number_of_Parameters;
		std::vector<std::string> parameter_names;
	};
#endif
