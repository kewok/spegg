#ifndef PENGUINS_H
#define PENGUINS_H

#include <species/inds_stochastic.h>

#include "Penguin_Parents.h"
#include <species/add_kids/neonates_class.h>

class Penguins : public inds_stochastic
	{
	public:
		Penguins(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
		void addKids();

		using inds_stochastic::update;
		void update(inds_stochastic **species);

	protected:
		void initialize_demes();
		void setPhenotype(int index, int n);
		void assignSex(int index, int n);
	};
#endif
