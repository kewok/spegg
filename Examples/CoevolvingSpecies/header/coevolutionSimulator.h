#ifndef COEVOLUTION_SIMULATOR_H
#define COEVOLUTION_SIMULATOR_H

#include <math/statistics_class.h>
#include <math/demographic_statistics_class.h>

#include <Simulation_Class.h>
#include "coevolvingSpecie.h"
#include "myConsumer.h"
#include "myResource.h"

class Coevolution_Simulator : public Simulation
	{
	public:
		Coevolution_Simulator();
		~Coevolution_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_predator_attack_ability;
		Statistics *stats_prey_defense_ability;
		Statistics *stats_prey_competitive_ability;
		DemographicStatistics *demographics_predator;
		DemographicStatistics *demographics_prey;

		void initialize_classes();
		static coevolvingSpecie *initialize_species(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);

		int nspecies;
	};

#endif
