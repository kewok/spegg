#ifndef MY_RESOURCE_H
#define MY_RESOURCE_H

#include "coevolvingSpecie.h"

class myResource : public coevolvingSpecie
	{
	public:
		myResource(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
	};

#endif
