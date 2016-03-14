#ifndef COEVOLVING_SPECIE_H
#define COEVOLVING_SPECIE_H

#include <species/inds_stochastic_migratory.h>
#include <species/add_kids/neonates_class.h>
#include "coevolvingSpeciesParents.h"
#include "deme_specific_data_class_coevolvingSpecies.h"

class coevolvingSpecie : public inds_stochastic_migratory
	{
	public:
		friend class update_coevolvingSpecies;
		friend class migrate_coevolvingSpecies;

		using inds_stochastic_migratory::update;
		using inds_stochastic_migratory::migrate;

		coevolvingSpecie(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);

		void addKids();
		void update(inds_stochastic **species);
		void migrate();

	protected:
		DemeSettings_coevolvingSpecies *demeParameters;
		void get_deme_data_vector(const char *parameter_name);
		void initialize_demes();
		void setPhenotype(int index, int n);
		void assignSex(int index, int n);		
	};

#endif 
