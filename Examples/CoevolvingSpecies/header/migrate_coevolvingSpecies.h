#ifndef MIGRATE_COEVOLVING_SPECIES_H
#define MIGRATE_COEVOLVING_SPECIES_H

#include <species/movement/migration_behavior.h>
#include "coevolvingSpecie.h"

#include <thrust/sequence.h>
#include <thrust/transform.h>

class migrate_coevolvingSpecies : public MigrationBehavior
	{
	coevolvingSpecie *species;
	// Constructor
	public:
		migrate_coevolvingSpecies(inds_stochastic_migratory *species) 
		 	{
			this->species = (coevolvingSpecie *) species;
			
			// Copy the constants 
			this->size = species->size;
			this->Number_of_Demes = species->Num_Demes;

			migration_offsets.resize(Number_of_Demes);
			thrust::sequence(migration_offsets.begin(), migration_offsets.end(), Number_of_Demes-1, Number_of_Demes);

			will_migrate.resize(size);
			thrust::fill(will_migrate.begin(), will_migrate.begin() + size, 0);

			migrant_destinations.resize(size);

			random_gen = species->gen;
			}

		void migrate();
	};

#endif
