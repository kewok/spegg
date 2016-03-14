#include <species/movement/migration_behavior.h>

void MigrationBehavior::determine_if_individuals_migrate(inds_stochastic_migratory *species)
	{
	/*
	* A generic template for the code that decides whether individuals will migrate during each time step. 
	*/
	thrust::host_vector<float> rand(size);
	for (int i=0; i < rand.size(); i++)
		{
		rand[i] = gsl_rng_uniform(species->gen);
		if (rand[i] < species->individual_migration_rate[i])
			{
			will_migrate[i] = 1;
			}
		}
	}

void MigrationBehavior::move_individuals(inds_stochastic_migratory *species)
	{
	thrust::host_vector<float> rand(size);

	for (int i=0; i < rand.size(); i++)
		{
		rand[i] = gsl_rng_uniform(species->gen);
		if (rand[i] < species->individual_migration_rate[i])
			{
			int inds_deme = species->deme[i];
			species->Migration_Matrix[inds_deme]->draw(species->gen, species->deme[i]);
			}
		}
	}
