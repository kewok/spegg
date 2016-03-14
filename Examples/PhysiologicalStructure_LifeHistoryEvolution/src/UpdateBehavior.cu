#include <species/update/updatebehavior.h>
#include "update_Fish.h"

UpdateBehavior *UpdateBehavior::create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID)
	{
	if (species_ID==0)
		return new update_Fish(species[0], habitat);
	}
