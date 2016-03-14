#include "update_Penguins.h"

UpdateBehavior * UpdateBehavior::create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID)
    {
    if (species_ID==0)
	   return new update_Penguins(species[0]);
    }
