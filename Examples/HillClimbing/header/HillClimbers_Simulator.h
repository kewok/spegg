#ifndef HILLCLIMBERS_SIMULATOR_H
#define HILLCLIMBERS_SIMULATOR_H

#include <Simulation_Class.h>
#include "HillClimbers.h"
#include <math/statistics_class.h>
#include <math/demographic_statistics_class.h>

class HillClimbers_Simulator : public Simulation
	{
	public:
		HillClimbers_Simulator();
		~HillClimbers_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_hillclimber_fecundity;
		Statistics *stats_hillclimber_mortality;
		DemographicStatistics *demographics;
		void initialize_classes();

		int nspecies;
	};
#endif
