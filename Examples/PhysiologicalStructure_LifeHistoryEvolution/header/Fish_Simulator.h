#ifndef FISH_SIMULATOR_H
#define FISH_SIMULATOR_H

#include <Simulation_Class.h>
#include "Fish.h"
#include "Fish_Habitat.h"
#include <math/statistics_class.h>
#include <math/demographic_statistics_class.h>

class Fish_Simulator : public Simulation
	{
	public:
		Fish_Simulator();
		~Fish_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_eggsize;
		Statistics *stats_fecundity;
		DemographicStatistics *demographics;
		Fish_Habitat *habitat;
		std::ofstream preyfile;

		void initialize_classes();

		int nspecies;
	};

#endif
