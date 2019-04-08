#include <species/deme_specific_data_class.h>

DemeSettings::DemeSettings(const char *filename, int species_ID)
	{
	read_in_parameters(filename, species_ID);
	specify_parameter_index();

	GeneticArchitecture = new DemeGeneticsSettings(filename, species_ID);
	}

void DemeSettings::specify_parameter_index()
	{
	for (int i=0; i < parameter_names.size(); i++)
		parameter_index[parameter_names[i]] = i;
	}

void DemeSettings::read_in_parameters(const char *filename, int species_ID)
	{
	/* 
	* A series of steps for reading the configuration file. Potentially this too can be broken up into different methods, with a single method for a single libconfig setting.
	*/
	Config cfg;

	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No deme config file." << std::endl;
		}
	catch(const ParseException &pex)
		{
		std::cerr << "Your " << pex.getFile() << " file is incorrectly specified. Make sure you check on or about line: " << pex.getLine() << " - " << pex.getError() << std::endl;
		}
	try
		{
		Number_of_Demes = cfg.lookup("number_of_demes");
		}
 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'number_of_demes' setting in configuration file." << std::endl;
		}

	const Setting& root = cfg.getRoot();

	// Read in the input

	const Setting &species_specification = root["species_data"];

	int nspecies = species_specification.getLength();
	
	species_specification[species_ID].lookupValue("number_of_parameters", Number_of_Parameters);
	
	parameter_names.resize(Number_of_Parameters);

	if (parameter_names.size() > 0)
		{
		if (species_specification[species_ID]["parameter_names"].getLength() != Number_of_Parameters)
			std::cout << "Length of parameter names differs from number of parameters." << std::endl;

		for (int i=0; i < species_specification[species_ID]["parameter_names"].getLength(); i++)
			{
			parameter_names[i] = species_specification[species_ID]["parameter_names"][i].c_str();
			}
		}
	else
		{
		std::cout << "warning: your species has no parameters in the deme configuration file." << std::endl;
		}

	// Read in the species_specific indices/values that apply across all demes
	try
		{
		species_specification[species_ID].lookupValue("number_of_species_specific_values", Number_of_Species_Specific_Values);
		}

 	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'number_of_species_specific_values' setting in configuration file." << std::endl;
		}

	species_specific_values_names.resize( Number_of_Species_Specific_Values );

	if (species_specific_values_names.size() > 0)
		{
		if (species_specification[species_ID]["species_specific_values_names"].getLength() != Number_of_Species_Specific_Values)
			std::cout << "number_of_species_specific_values doesn't match length of species specific values list in config file" << std::endl;

		for (int i=0; i < species_specification[species_ID]["species_specific_values_names"].getLength(); i++)
			{
			species_specific_values_names[i] = species_specification[species_ID]["species_specific_values_names"][i].c_str();
			}
		}
	else
		{
		std::cout << "warning: your species has no species-specific values in the deme configuration file." << std::endl;
		}


	for (int i=0; i < Number_of_Species_Specific_Values; i++)
		{
		float val = 0;
		species_specification[species_ID].lookupValue(species_specific_values_names[i], val);
		species_specific_values[species_specific_values_names[i]] = val;
		}

	if (species_specification[species_ID].exists("demes_specifications"))
		{
		const Setting &deme_specifications = species_specification[species_ID]["demes_specifications"];
	
		deme_wide_parameters = new thrust::device_vector<float>[Number_of_Parameters];

		for (int i=0; i < Number_of_Parameters; i++)
			deme_wide_parameters[i].resize(Number_of_Demes);

		for (int i=0; i < parameter_names.size(); i++)
			{
			for (int j=0; j < Number_of_Demes; j++)
				{
				const Setting &deme_values = deme_specifications[j];
				float val = 0;
				deme_values.lookupValue(parameter_names[i], val);
				deme_wide_parameters[i][j] = val;
				}
			}
		}
	else
		{
		std::cerr << "Warning: No 'deme_specifications' setting in configuration file." << std::endl;
		}
	}

thrust::device_ptr<float> DemeSettings::get_vector_ptr(const char *parameter_name)
	{
	return(&deme_wide_parameters[parameter_index[parameter_name]][0]);
	}

int DemeSettings::check_number_of_demes()
	{
	return(Number_of_Demes);
	}


bool DemeSettings::does_parameter_exist(const char *parameter_name)
	{
	if(std::find(parameter_names.begin(), parameter_names.end(), parameter_name) != parameter_names.end())
		{
		return true;
		}
		
	else
		{
		return false;
		}
	}

