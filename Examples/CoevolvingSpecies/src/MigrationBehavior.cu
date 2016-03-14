#include <species/movement/migration_behavior.h>
#include "migrate_coevolvingSpecies.h"

MigrationBehavior *MigrationBehavior::create_migrationBehavior(inds_stochastic_migratory *species)
	{
	return new migrate_coevolvingSpecies(species);
	}

