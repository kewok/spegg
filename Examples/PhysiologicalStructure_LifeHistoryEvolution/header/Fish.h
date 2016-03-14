#ifndef FISH_H
#define FISH_H

#include <species/inds_stochastic.h>
#include "Fish_Habitat.h"
#include "Fish_Parents.h"
#include <species/add_kids/neonates_class.h>

class Fish : public inds_stochastic
	{
	public:
		using inds_stochastic::update;
		Fish(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
		void addKids();
		void update(inds_stochastic **species, environment *habitat);

	protected:
		void get_deme_data_vector(const char *parameter_name);

		void initialize_demes();

		void setPhenotype(int index, int n);
		void assignSex(int index, int n);
	};

#endif

