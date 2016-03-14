#ifndef PENGUIN_SIMULATOR_H
#define PENGUIN_SIMULATOR_H

#include <Simulation_Class.h>
#include "Penguins.h"
#include <math/statistics_class.h>

class Penguin_Drift_Simulator : public Simulation
	{
	public:
		Penguin_Drift_Simulator();
		~Penguin_Drift_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_penguins;
		void initialize_classes();

		int nspecies;
	};
#endif
