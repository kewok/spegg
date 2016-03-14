#include "deme_specific_data_class_coevolvingSpecies.h"

DemeSettings_coevolvingSpecies ::DemeSettings_coevolvingSpecies (const char *filename, int species_ID) : DemeSettings(filename, species_ID)
	{
	read_in_interactions(filename, species_ID);
	}

void DemeSettings_coevolvingSpecies ::read_in_interactions(const char *filename, int species_ID)
	{
	Config cfg;

	try
		{
		cfg.readFile(filename);
		}
	catch(const FileIOException &fioex)
		{
		std::cerr << "No deme config file." << std::endl;
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

	/* read in the interactions */
	const Setting &interacting_species_indices = species_specification[species_ID]["interacting_species"];

	interacting_species.resize(interacting_species_indices.getLength());

	for (int i=0; i < species_specification[species_ID]["interacting_species"].getLength(); i++)
		{
		int val = 0;
		val = interacting_species_indices[i];
		interacting_species[i] = val;
		}

	const Setting &interaction_phenotype_values = species_specification[species_ID]["interaction_phenotype_indices"];

	for (int i=0; i < interacting_species_indices.getLength(); i++)
		{
		interaction_phenotype_indices[interacting_species[i]] = interaction_phenotype_values[i];
		}

	deme_wise_interaction_effects_on_fecundity = new thrust::host_vector<float>[interacting_species_indices.getLength()];
	deme_wise_interaction_effects_on_survivorship = new thrust::host_vector<float>[interacting_species_indices.getLength()];

	const Setting &deme_specifications = species_specification[species_ID]["demes_specifications"];

	int Number_of_Demes = deme_specifications.getLength();

	for (int i=0; i < interacting_species_indices.getLength(); i++)
		{
		deme_wise_interaction_effects_on_fecundity[i].resize(Number_of_Demes);
		deme_wise_interaction_effects_on_survivorship[i].resize(Number_of_Demes);
		}

	for (int i=0; i < interacting_species_indices.getLength(); i++)
		{
		for (int j=0; j < Number_of_Demes; j++)
			{
			const Setting &deme_wise_interactions = deme_specifications[j]["interaction_effects_on_fecundity"];
			deme_wise_interaction_effects_on_fecundity[i][j] = deme_wise_interactions[i];
			}

		for (int j=0; j < Number_of_Demes; j++)
			{
			const Setting &deme_wise_interactions = deme_specifications[j]["interaction_effects_on_survivorship"];
			deme_wise_interaction_effects_on_survivorship[i][j] = deme_wise_interactions[i];
			}
		}
	}

