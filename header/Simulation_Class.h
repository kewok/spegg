#ifndef SIM_H
#define SIM_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include <libconfig.h++>

class Simulation
	{
	public: 
		Simulation();
		~Simulation();
		virtual void run() = 0;
		
	protected:
		// Input data
		int threadID;
		int nsteps;
		int demes;
		int initpop;
		int maxpop;
		int intra_step_time_steps;
		int nloci;
		int nphenotypes;
		int seed;

		int num_biotic_variables;
		int num_abiotic_variables;
		std::vector<int> steps_to_output_individuals_csv;

		// methods
		void read_simulation_settings();
		virtual void initialize_classes() = 0;
	};

#endif
