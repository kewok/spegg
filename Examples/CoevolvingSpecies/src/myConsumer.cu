#include "myConsumer.h"

myConsumer::myConsumer(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : coevolvingSpecie(size_val, maxsize_val, seed_val, ndemes, 0)
	{
	initialize_demes();

	for (int i=0; i < nloci; i++)
		{		
		draw_gaussian(size, 0.0, 1.0, fgenotype[i], gen);
		draw_gaussian(size, 0.0, 1.0, mgenotype[i], gen);

		// Just make the genotypes positive 
		for (int j=0; j < size; j++)
			{
			fgenotype[i][j] = fabs(fgenotype[i][j]);
			mgenotype[i][j] = fabs(mgenotype[i][j]);
			}
		}
	

	// To start, assign odd numbered individuals to be male, even numbered individuals to be female
	
	thrust::host_vector<int> twos(size);
	thrust::fill(twos.begin(), twos.end(), 2);
	thrust::transform(id.begin(), id.begin() + size, twos.begin(), sex.begin(), thrust::modulus<int>());
	
	//Set phenotype
	setPhenotype(0, size);	
	}


