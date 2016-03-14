#include <species/update/updatebehavior.h>
#include "update_coevolvingSpecies.h"

UpdateBehavior *UpdateBehavior::create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID)
	{
	return new update_coevolvingSpecies(species, species_ID);
	}
