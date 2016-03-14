#include <Simulation_Class.h>

Simulation::Simulation()
	{
	read_simulation_settings();
	}
	
void Simulation::read_simulation_settings()
	{
	ConfigFile config("Simulation.conf");
	nsteps = config.read<int>( "n_timesteps" );
	demes = config.read<int>( "ndemes" );
	intra_step_time_steps = config.read<int>( "intra_step_time_steps" );

	seed = config.read<int>("random_seed");

	num_biotic_variables = config.read<int>( "num_biotic_variables" );
	num_abiotic_variables = config.read<int>( "num_abiotic_variables" );
	}

Simulation::~Simulation()
	{
	}
