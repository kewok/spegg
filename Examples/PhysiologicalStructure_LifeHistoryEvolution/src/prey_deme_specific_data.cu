#include "prey_deme_specific_data.h"

PreyDemeSpecificData::PreyDemeSpecificData(const char *filename, int prey_index)
	{
	read_in_prey_parameters(filename, prey_index);
	specify_parameter_index();
	}

void PreyDemeSpecificData::specify_parameter_index()
	{
	for (int i=0; i < parameter_names.size(); i++)
		parameter_index[parameter_names[i]] = i;
	}

void PreyDemeSpecificData::read_in_prey_parameters(const char *filename, int prey_index)
	{
	Config cfg;
	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No prey config file." << std::endl;
		}
	
	try
		{
		Number_of_Parameters = cfg.lookup("number_of_parameters");
		}
	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "No 'number_of_parameters' setting in prey configuration file." << std::endl;
		}

	try
		{
		const Setting& root = cfg.getRoot();
		parameter_names.resize(Number_of_Parameters);

		if (root["parameter_names"].getLength() != Number_of_Parameters)
			std::cout << "Length of parameter names differs from number of parameters\n" << std::endl;

		for (int i=0; i < root["parameter_names"].getLength(); i++)
				{
				parameter_names[i] = root["parameter_names"][i].c_str();
				}

		const Setting &deme_specifications = root["prey_initialization"][prey_index];

		Number_of_Demes = deme_specifications.getLength();

		deme_wide_parameters = new thrust::device_vector<float>[Number_of_Parameters];

		for (int i=0; i < Number_of_Parameters; i++)
			deme_wide_parameters[i].resize(Number_of_Demes);

		for (int i=0; i < root["parameter_names"].getLength(); i++)
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
	
	catch(const SettingNotFoundException &nfex)
		{
		std::cerr << "Your prey configuration file does not work." << std::endl;
		}
	}

thrust::device_ptr<float> PreyDemeSpecificData::get_vector_ptr(const char *parameter_name)
	{
	return(&deme_wide_parameters[parameter_index[parameter_name]][0]);
	}
