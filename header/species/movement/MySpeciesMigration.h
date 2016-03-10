// An example template for how a migrate_mySpecies function might work.

#include "migration_behavior.h"
// Include your MySpecies class header file here

class MySpeciesMigration : public MigrationBehavior
{
	class MySpecies *species;
	class MyHabitat **habitat;

	// constructor; in this example, environmentally cued dispersal is potentially allowed by including MyHabitat **habitat in the arguments.
	public:
	MySpeciesMigration(MySpecies *species, MyHabitat **habitat)
		{
		this->species = species;
		this->habitat = habitat;
		this->size = species->size;
		this->Number_of_Demes = species->Num_Demes;
		this->random_gen = species->gen;

		will_migrate.resize(size);
		thrust::fill(will_migrate.begin(), will_migrate.begin() + size, 0);

		migrant_destinations.resize(size);

		// The rationale for the migration_offsets vector is obscure. It has to do with subpop_thrust_prob_table whose code is a mess. It's also not clear it should be a member of MySpeciesMigration. It should in fact probably be a member of subpop_thrust_prob_table.
		migration_offsets.resize(Number_of_Demes);
		thrust::sequence(migration_offsets.begin(), migration_offsets.end(), Number_of_Demes-1, Number_of_Demes);
		}

	// busines	
	void simulate_migration();
}
