#ifndef FISH_HABITAT_H
#define FISH_HABITAT_H

#include <species/inds.h>
#include <environ/environment.h>
#include "prey_class.h"

class Fish_Habitat : public environment
	{
	public:
		Fish_Habitat(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes, int num_time_steps);
		void update();

		prey **prey_array; // An array of prey classes
	
		int intra_annual_time_steps;
		int current_time_step;
	protected:
		void initialize_abiotic_variables(const char *filename);
	};

#endif
