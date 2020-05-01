#include "Simulation_Class.h"

using namespace libconfig;

Simulation::Simulation()
	{
	read_simulation_settings();
	}
	
void Simulation::read_simulation_settings()
	{
	Config cfg;

	try
		{
		cfg.readFile("Simulation.conf");
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No Simulation.conf file." << std::endl;
		}
	catch(const ParseException &pex)
		{
		std::cerr << "Your " << pex.getFile() << " file is incorrectly specified. Make sure you check on or about line: " << pex.getLine() << " - " << pex.getError() << std::endl;
		}
	try
		{
		nsteps = cfg.lookup("n_timesteps");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'n_timesteps' setting in configuration file." << std::endl;
		}
	try
		{
		demes = cfg.lookup("ndemes");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'ndemes' setting in configuration file." << std::endl;
		}
	try
		{
		num_biotic_variables = cfg.lookup("num_biotic_variables");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'num_biotic_variables' setting in configuration file." << std::endl;
		}
	try
		{
		num_biotic_variables = cfg.lookup("num_abiotic_variables");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'num_abiotic_variables' setting in configuration file." << std::endl;
		}
	// Optional entries
	try
		{
		intra_step_time_steps = cfg.lookup("intra_step_time_steps");
		}
	catch(const SettingNotFoundException &nfex)
		{
		// Ignore.
		}
	try
		{
		seed = cfg.lookup("random_seed");
		}
	catch(const SettingNotFoundException &nfex)
		{
		// Ignore.
		}
	// Individual csv recording times
	try
		{
		const Setting& root = cfg.getRoot();
		const Setting &output_csv_steps = root["output_csv_steps"];
		for (int i=0; i < output_csv_steps.getLength(); i++)
			{
			steps_to_output_individuals_csv.push_back(output_csv_steps[i]);
			std::cout << steps_to_output_individuals_csv[i] << std::endl;
			}
		}
	catch(const SettingNotFoundException &nfex)
		{
		// Ignore.
		}
	}

Simulation::~Simulation()
	{
	}
