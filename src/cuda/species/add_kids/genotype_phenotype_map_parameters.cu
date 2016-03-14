#include <species/add_kids/genotype_phenotype_map_parameters.h>

// For simplicity, right now assume that phenotypes are arranged correctly in the deme_config.txt (this might in any event be the responsibility of the config file generating program, not something these classes should have to worry about. 

GenotypePhenotypeMapParameters::GenotypePhenotypeMapParameters(const char *filename, int species_ID, int phenotype_index, std::vector<std::string> &loci_names)
	{
	this->phenotype_index = phenotype_index;
	read_in_data(filename, species_ID);
	specify_parameter_index();
	}


void GenotypePhenotypeMapParameters::specify_parameter_index()
	{
	for (int i=0; i < Names_of_Genotype_Phenotype_Map_Parameters.size(); i++)
		parameter_index[Names_of_Genotype_Phenotype_Map_Parameters[i]] = i;
	}

void GenotypePhenotypeMapParameters::read_in_data(const char *filename, int species_ID)
	{
	Config cfg;
	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No genetic config file. \nyou failure. much dishonor to family name. suggest immediate death by seppuku." << std::endl;
		}
	
	const Setting& root = cfg.getRoot();

	const Setting &species_genetic_specification = root["species_genetics"];
	
	const Setting &phenotype_specification = species_genetic_specification[species_ID]["phenotype_specifications"][phenotype_index];

	phenotype_specification.lookupValue("number_of_genotype_phenotype_map_parameters", Number_of_Parameters);

	if (Number_of_Parameters > 0)
		{
		Names_of_Genotype_Phenotype_Map_Parameters.resize(Number_of_Parameters);

		if (phenotype_specification["names_of_genotype_phenotype_map_parameters"].getLength() != Number_of_Parameters)
			std::cout << "Length of parameter names differs from number of parameters\n" << std::endl;

		for (int i=0; i < phenotype_specification["names_of_genotype_phenotype_map_parameters"].getLength(); i++)
			{
			Names_of_Genotype_Phenotype_Map_Parameters[i] = phenotype_specification["names_of_genotype_phenotype_map_parameters"][i].c_str();
			}

		const Setting &parameter_values = phenotype_specification["genotype_phenotype_map_parameters"];

		Number_of_Demes = parameter_values.getLength();

		deme_specific_parameters = new thrust::host_vector<float>[Number_of_Parameters];

		for (int i=0; i < Number_of_Parameters; i++)
			deme_specific_parameters[i].resize(Number_of_Demes);

		for (int i=0; i < Names_of_Genotype_Phenotype_Map_Parameters.size(); i++)
			{
			for (int j=0; j < Number_of_Demes; j++)
				{
				const Setting &deme_values = parameter_values[j]; // go to deme j
				float val = 0;
				deme_values.lookupValue(Names_of_Genotype_Phenotype_Map_Parameters[i], val);
				deme_specific_parameters[i][j] = val;
				}
			}
		}
	}

float *GenotypePhenotypeMapParameters::get_vector_ptr(const char *parameter_name)
	{
	return(&deme_specific_parameters[parameter_index[parameter_name]][0]);
	}


