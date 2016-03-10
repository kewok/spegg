#ifndef UPDATE_BEHAVIOR_H
#define UPDATE_BEHAVIOR_H

#include "inds_stochastic.h"
#include "survivorship_kernel_functors.h"
#include "environment.h"
#include <thrust/sequence.h>

//virtual interface; 

class UpdateBehavior
	{
	public:
		/* the factory */
		static UpdateBehavior *create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID);

		/* the actual updating */
		virtual void update()=0;

		virtual ~UpdateBehavior() {};

		/* Generic functionality for simulating mortality by changing the vital state variable according to a bernoulli RV */
		void determine_mortality(inds_stochastic *species);
	};

#endif
