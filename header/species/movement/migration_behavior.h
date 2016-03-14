#ifndef MIGRATION_BEHAVIOR_H
#define MIGRATION_BEHAVIOR_H

#include <species/inds_stochastic_migratory.h>
#include <math/thrust_prob_table_demes.h>
#include <species/movement/migration_kernel_functors.h>

class MigrationBehavior
	{
	public:
		/* the factory */
		static MigrationBehavior *create_migrationBehavior(inds_stochastic_migratory *species);
		
		/* the actual migration */
		virtual void migrate()=0;

		virtual ~MigrationBehavior() {};

		/* convenient constants used throughout */
		int size, Number_of_Demes;
		curandGenerator_t random_gen;

		/* helper vectors */
		thrust::device_vector<int> migration_offsets;
		thrust::device_vector<int> migrant_destinations;
		thrust::device_vector<float> will_migrate;

		/* main migration functionality. move_individuals is the primary interface */
		void move_individuals(inds_stochastic_migratory *species);
		void determine_destination(inds_stochastic_migratory *species);
		void determine_if_individuals_migrate(inds_stochastic_migratory *species);
	};

#endif
