#include "myConsumer.h"

myConsumer::myConsumer(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : coevolvingSpecie(size_val, maxsize_val, seed_val, ndemes, 0)
	{
	initialize_demes();

	for (int i=0; i < nloci; i++)
		{		
		draw_gaussian(size, 0.0, 1.0, fgenotype[i], gen);
		draw_gaussian(size, 0.0, 1.0, mgenotype[i], gen);

		// Just make the genotypes positive 
		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(fgenotype[i].begin())),
				 thrust::make_zip_iterator(thrust::make_tuple(fgenotype[i].begin() + size)), 
				 fabs_functor());
		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(mgenotype[i].begin())),
				 thrust::make_zip_iterator(thrust::make_tuple(mgenotype[i].begin() + size)), 
				 fabs_functor());
		}
	

	// To start, assign odd numbered individuals to be male, even numbered individuals to be female
	
	thrust::device_vector<int> twos(size);
	thrust::fill(twos.begin(), twos.end(), 2);
	thrust::transform(id.begin(), id.begin() + size, twos.begin(), sex.begin(), thrust::modulus<int>());
	
	//Set phenotype
	setPhenotype(0, size);	
	}


