#ifndef MY_CONSUMER_H
#define MY_CONSUMER_H

#include "coevolvingSpecie.h"

class myConsumer : public coevolvingSpecie
	{
	public:
		myConsumer(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
	};

#endif
