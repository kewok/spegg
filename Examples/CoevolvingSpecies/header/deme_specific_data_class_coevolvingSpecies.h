#ifndef DEME_SPECIFIC_DATA_CLASS_COEVOLVING_SPECIES_H
#define DEME_SPECIFIC_DATA_CLASS_COEVOLVING_SPECIES_H

// Requires libconfig++; after installing, add /usr/local/lib to path via LD_LIBRARY_PATH as well, followed by sudo ldconfig!

#include <species/deme_specific_data_class.h>

using namespace libconfig;

class DemeSettings_coevolvingSpecies : public DemeSettings
{
public: 
	DemeSettings_coevolvingSpecies(const char *filename, int species_ID);

	thrust::host_vector<int> interacting_species; // species indices with which this species interacts

	// The map interaction_phenotype_indices is designed to take as its argument the species index with which the interaction is occurring, and return the index of the phenotype corresponding governing this interaction.
	std::map<int, int> interaction_phenotype_indices;

	thrust::host_vector<float> *deme_wise_interaction_effects_on_fecundity;
	thrust::host_vector<float> *deme_wise_interaction_effects_on_survivorship;

	void read_in_interactions(const char *filename, int species_ID);
};

#endif
