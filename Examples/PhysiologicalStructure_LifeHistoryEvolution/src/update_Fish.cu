#include "update_Fish.h"
#include "Fish_Habitat.h"

void update_Fish::update()
	{
	float *resource_1 = raw_pointer_cast(&(habitat->biotic_variables[0][0]));
	float *resource_2 = raw_pointer_cast(&(habitat->biotic_variables[1][0]));

	grand_consumption_functor eat_N_grow(consumption_allometric_scalar, gamma, ontogenetic_niche_shift_constant, ontogenetic_niche_shift_coefficient, length_weight_conversion_coefficient, length_weight_conversion_exponent, resource_1, 	resource_1_maximum, resource_2, resource_2_maximum, handling_time_resource1, handling_time_resource2, functional_response_scalar_resource1, functional_response_scalar_resource2, alpha_g, mature_maximum_condition, juvenile_maximum_condition, M_sizes_at_maturity, F_sizes_at_maturity);

	thrust::device_vector<float> resource_1_consumed(size);
	thrust::device_vector<float> resource_2_consumed(size);

	thrust::device_vector<float> rand(size);
	float *rand_ptr = raw_pointer_cast(&rand[0]);
		
	thrust::device_vector<float> temp_eaten(Number_of_Demes);
	thrust::device_vector<int> temp_demes(Number_of_Demes);
	
	//Set up mortality functor
	calculate_mortality mortality_functor(effect_of_starvation_ptr, size_dependent_mortality_constant_ptr, size_dependent_mortality_coefficient_ptr);

	for (int Time_Step=0; Time_Step < intra_annual_time_steps; Time_Step++)
		{
		// Draw the random numbers for somatic growth
		curandGenerateUniform(species->gen, rand_ptr, size);

		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin(), species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin(), species->deme.begin(), species->phenotype[RESOURCE_LIMITATION_PHENOTYPE].begin(), species->sex.begin(), species->status.begin(), resource_1_consumed.begin(), resource_2_consumed.begin(), rand.begin())),
				 thrust::make_zip_iterator(thrust::make_tuple(species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size, species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin() + size, species->deme.begin() + size, species->phenotype[RESOURCE_LIMITATION_PHENOTYPE].begin() + size,  species->sex.begin() + size, species->status.begin() + size, resource_1_consumed.begin() + size, resource_2_consumed.begin() + size, rand.end() + size)), 
				eat_N_grow);
	

		// Draw the random numbers for survivorship
		curandGenerateUniform(species->gen, rand_ptr, size);

		// Apply survivorship
		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin(), species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin(), species->phenotype[MORTALITY_PHENOTYPE].begin(), species->deme.begin(), species->age.begin(), species->status.begin(), rand.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size, species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin() + size, species->phenotype[MORTALITY_PHENOTYPE].begin() + size, species->deme.begin() + size,  species->age.begin() + size, species->status.begin() + size, rand.begin() + size)),
				mortality_functor);
			

		// Calculate the values for the environmental variables
		cudaDeviceSynchronize();
		reduce_by_key(species->deme.begin(), species->deme.begin() + size, resource_1_consumed.begin(), temp_demes.begin(), temp_eaten.begin());
		thrust::scatter(temp_eaten.begin(), temp_eaten.end(), temp_demes.begin(), habitat->effect_of_inds_on_biotic_variable[0].begin());
		reduce_by_key(species->deme.begin(), species->deme.begin() + size, resource_2_consumed.begin(), temp_demes.begin(), temp_eaten.begin());
		thrust::scatter(temp_eaten.begin(), temp_eaten.end(), temp_demes.begin(), habitat->effect_of_inds_on_biotic_variable[1].begin());

		habitat->update();
		}

	update_vital_rates();
	}

void update_Fish::prepare_survivorship_constants_pointers()
	{
	effect_of_starvation_ptr = raw_pointer_cast(species->demeParameters->get_vector_ptr("effect_of_starvation"));
	size_dependent_mortality_constant_ptr = raw_pointer_cast(species->demeParameters->get_vector_ptr("size_dependent_mortality_constant"));
	size_dependent_mortality_coefficient_ptr = raw_pointer_cast(species->demeParameters->get_vector_ptr("size_dependent_mortality_coefficient"));
	}

void update_Fish::prepare_growth_constants_pointers()
	{
	consumption_allometric_scalar = raw_pointer_cast(species->demeParameters->get_vector_ptr("conversion_efficiency"));
	gamma = raw_pointer_cast(species->demeParameters->get_vector_ptr("gamma_exponent"));
	ontogenetic_niche_shift_constant = raw_pointer_cast(species->demeParameters->get_vector_ptr("ontogenetic_niche_shift_constant"));
	ontogenetic_niche_shift_coefficient = raw_pointer_cast(species->demeParameters->get_vector_ptr("ontogenetic_niche_shift_coefficient"));
	length_weight_conversion_coefficient = raw_pointer_cast(species->demeParameters->get_vector_ptr("length_weight_conversion_coefficient"));
	length_weight_conversion_exponent = raw_pointer_cast(species->demeParameters->get_vector_ptr("length_weight_conversion_exponent"));

	resource_1_maximum = raw_pointer_cast(&(habitat->prey_array[0]->prey_maximum_abundance[0]));
	resource_2_maximum = raw_pointer_cast(&(habitat->prey_array[1]->prey_maximum_abundance[0]));
	handling_time_resource1 = raw_pointer_cast(species->demeParameters->get_vector_ptr("handling_time"));
	handling_time_resource2 = raw_pointer_cast(species->demeParameters->get_vector_ptr("handling_time"));

	functional_response_scalar_resource1 = raw_pointer_cast(species->demeParameters->get_vector_ptr("functional_response_numerator"));
	functional_response_scalar_resource2 = raw_pointer_cast(species->demeParameters->get_vector_ptr("functional_response_numerator"));

	/* constants governing growth rates */
	alpha_g = raw_pointer_cast(species->demeParameters->get_vector_ptr("alpha_g"));
	mature_maximum_condition = raw_pointer_cast(species->demeParameters->get_vector_ptr("mature_maximum_condition"));
	juvenile_maximum_condition = raw_pointer_cast(species->demeParameters->get_vector_ptr("juvenile_maximum_condition"));

	M_sizes_at_maturity = raw_pointer_cast(species->demeParameters->get_vector_ptr("M_sizes_at_maturity"));
	F_sizes_at_maturity = raw_pointer_cast(species->demeParameters->get_vector_ptr("F_sizes_at_maturity"));
	}

void update_Fish::update_vital_rates()
	{
	update_fecundity();
	}

void update_Fish::update_fecundity()
	{
	//Set up reproductive contribution functor.
	float *juvenile_maximum_conditions_ptr = raw_pointer_cast(species->demeParameters->get_vector_ptr("juvenile_maximum_condition"));

	update_female_fecundity_functor fecundity(juvenile_maximum_conditions_ptr);

	thrust::device_vector<int> indices(size);
	thrust::sequence(indices.begin(), indices.begin() + size);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->sex.begin(), species->deme.begin(), species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin(), species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin(), species->phenotype[EGGSIZE_PHENOTYPE].begin(), species->phenotype[FECUNDITY_PHENOTYPE].begin(), indices.begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(species->sex.begin() + size, species->deme.begin() + size, species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE].begin() + size, species->phenotype[REVERSIBLE_MASS_PHENOTYPE].begin() + size, species->phenotype[EGGSIZE_PHENOTYPE].begin() + size, species->phenotype[FECUNDITY_PHENOTYPE].begin() + size, indices.end())),
			 fecundity);
	}
