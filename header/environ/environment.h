#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <util/thrust_functors.h>
#include <map>

#include <libconfig.h++>

using namespace libconfig;

class environment
	{
	public:
		environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes);
		void initialize_abiotic_variables(const char *filename);
		~environment();

		//Random number generator
		curandGenerator_t gen;

		//Misc data ints
		int seed;
		int ndemes;
		int nbiotic_vars;
		int nabiotic_vars;

		//void update(int num_biotic_variables, int month);
		void record_habitat_details(const char *output_file_environment);
		thrust::device_ptr<float> get_abiotic_vector_ptr(const char *abiotic_variable_name);


		// Data vectors - effect_of_inds_on_biotic_variables[prey][deme]
		thrust::device_vector<float> *biotic_variables;
		thrust::device_vector<float> *effect_of_inds_on_biotic_variable;

		virtual void update() = 0;
		// Gradually move more of the other members into protected
	protected:
		std::vector<std::string> abiotic_variable_names;
		thrust::device_vector<float> *abiotic_variables;

		std::map<std::string, int> abiotic_variable_indices;
	};
#endif
