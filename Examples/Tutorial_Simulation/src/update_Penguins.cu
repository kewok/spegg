#include "update_Penguins.h"

void update_Penguins::update()
	{
	determine_mortality(species);

	// Change the crown color
	update_color();

	//Increment age for all survivors.
	thrust::transform_if(species->age.begin(), species->age.begin() + size, species->status.begin(), species->age.begin(), unary_plus<int>(1), thrust::identity<int>());
	}

void update_Penguins::update_color()
	{
	int CROWN_COLOR_INDEX = species->demeParameters->species_specific_values["CROWN_COLOR_INDEX"];

	// Step 1. Cast the thrust vectors into arrays 
	float *crown_colors = raw_pointer_cast(&species->phenotype[CROWN_COLOR_INDEX][0]);
	float *target_crown_colors = raw_pointer_cast(&species->demeParameters->get_vector_ptr("TARGET_CROWN_COLOR")[0]);
	
	// Step 2. Set up functor
	crown_color_updater crown_color_functor(crown_colors, target_crown_colors);

    	// Step 3. Create a vector which, for each individual, stores the crown color decay rate associated with their deme.
    		// Step 3a. Create a one-time vector, demewise_color_decays that copies the deme-specific color decay rates.
	thrust::device_vector<float> demewise_color_decays(species->Num_Demes);
	thrust::copy(species->demeParameters->get_vector_ptr("CROWN_COLOR_DECAY"), species->demeParameters->get_vector_ptr("CROWN_COLOR_DECAY") + species->Num_Demes, demewise_color_decays.begin());

    		// Step 3b. Create a vector, color_decays of length size (i.e., the number of living individuals) which stores, for each individual, the crown color decay rate associated with their deme.
	thrust::device_vector<float> color_decays(size);
	amplify_float(demewise_color_decays, species->deme_sizes, color_decays);

    	// Step 4. specify the indices of individuals
	thrust::device_vector<int> individuals(size);
	thrust::sequence(individuals.begin(), individuals.begin() + size, 0);

	// Step 5. Perform greying operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals.begin(), species->deme.begin(), color_decays.begin())),
	                 thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + size, species->deme.begin() + size, color_decays.begin() + size)),
	                 crown_color_functor);
	}
