#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <gsl/gsl_rng.h>
#include <thrust/host_vector.h>
#include <util/thrust_functors.h>
#include <map>

#include <libconfig.h++>

using namespace libconfig;

class environment
{
	public:
		environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes);
		~environment();

		//Random number generator
		gsl_rng* gen;

		//Misc data ints
		int seed;
		int ndemes;
		int nbiotic_vars;
		int nabiotic_vars;

		//void update(int num_biotic_variables, int month);
		void record_habitat_details(const char *output_file_environment);
		float *get_abiotic_vector_ptr(const char *abiotic_variable_name);


		// Data vectors - effect_of_inds_on_biotic_variables[prey][deme]
		thrust::host_vector<float> *biotic_variables;
		thrust::host_vector<float> *effect_of_inds_on_biotic_variable;

		virtual void update() = 0;
		// Gradually move more of the other members into protected
	protected:
		std::vector<std::string> abiotic_variable_names;
		thrust::host_vector<float> *abiotic_variables;

		std::map<std::string, int> abiotic_variable_indices;
};
#endif
