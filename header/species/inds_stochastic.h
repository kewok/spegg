#ifndef STOCHASTIC_INDS_H
#define STOCHASTIC_INDS_H

#include <species/inds.h>
#include <math/random_variables_functions.h>

class inds_stochastic : public inds
	{
	public:
		inds_stochastic(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
		curandGenerator_t gen;
	
		virtual void update(inds_stochastic **species) {};
		virtual void update(inds_stochastic **species, environment *habitat) {};
		virtual void update(inds_stochastic **species, environment *habitat, int intra_step_time_steps) {};
		virtual void update(inds_stochastic **species, environment *habitat, int intra_step_time_steps, int current_time_step) {};

		virtual void addKids() {};
		virtual void addKids(environment *habitat) {};

		int seed;

		int MORTALITY_PHENOTYPE_INDEX;
		int FECUNDITY_PHENOTYPE_INDEX;
	};

#endif
