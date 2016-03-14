#include "update_Fish.h"
#include "Fish_Habitat.h"

// Update fecundity according to (reversible-irreversible*max_cond)/(max_cond*eggsize + egg_size)

void update_female_fecundity_functor(float &maximum_condition, float &irreversible_mass, float &reversible_mass, float &egg_size, float &fecundity)
	{
	fecundity = 0; // make sure there is no carry-over from the previous time step.
	fecundity = (reversible_mass - irreversible_mass*maximum_condition)/(maximum_condition*egg_size + egg_size);

	// because we are dealing with the number of kids, convert to integers
	int tempVal = (int) fecundity;
	fecundity = (float) tempVal;

	float surplus_mass_used_in_reproduction = (reversible_mass - irreversible_mass*maximum_condition);

	/* make sure that you deduct the amount used from reproductive females */
	if (fecundity > 0)
		{
		reversible_mass -= surplus_mass_used_in_reproduction;
		}
	else
		fecundity = 0 ;
	};

void update_Fish::update()
	{
	thrust::host_vector<float> resource_1_consumed(species->size);
	thrust::host_vector<float> resource_2_consumed(species->size);
	
	for (int i=0; i < species->Num_Demes; i++)
		{
		habitat->effect_of_inds_on_biotic_variable[0][i] = 0;
		habitat->effect_of_inds_on_biotic_variable[1][i] = 0;
		}
		
	thrust::host_vector<float> temp_eaten(Number_of_Demes);
	thrust::host_vector<int> temp_demes(Number_of_Demes);
	
	// Define variables:
	float xt0, yt0, xt1, yt1, kappa, maturity_threshold,  available_resources_1,  available_resources_2, common_satiation_resource_1, common_satiation_resource_2, satiation_resource_1, satiation_resource_2, prey_eaten_1, prey_eaten_2, maximum_consumption, percent_adult_resource_eaten, percent_juvenile_resource, Eg, rand;

	float functional_response_scalar_resource1, functional_response_scalar_resource2, handling_time_resource1, handling_time_resource2, F_sizes_at_maturity, M_sizes_at_maturity, mature_maximum_condition, juvenile_maximum_condition, conversion_efficiency, ons_constant, ons_coefficient, lw_coefficient, lw_exponent, effect_of_starvation, size_dependent_mortality_coefficient, size_dependent_mortality_constant, gamma_val, alpha_val;
	for (int Time_Step=0; Time_Step < intra_annual_time_steps; Time_Step++)
		{
		int offset = 0;
		for (int j=0; j < species->Num_Demes; j++)
			{
			if (j > 0)
				{
				offset += species->deme_sizes[j-1];
				}

			int individuals_deme = j;
			functional_response_scalar_resource1 = *(functional_response_scalar_resource1_ptr + individuals_deme);
			functional_response_scalar_resource2 = *(functional_response_scalar_resource2_ptr + individuals_deme) ;		
			handling_time_resource1= *(handling_time_resource1_ptr  + individuals_deme);
			handling_time_resource2 = *(handling_time_resource2_ptr  + individuals_deme); 

			F_sizes_at_maturity = *(F_sizes_at_maturity_ptr  + individuals_deme); 
			M_sizes_at_maturity = *(M_sizes_at_maturity_ptr  + individuals_deme); 
			mature_maximum_condition = *(mature_maximum_condition_ptr  + individuals_deme); 
			juvenile_maximum_condition = *(juvenile_maximum_condition_ptr  + individuals_deme); 
			conversion_efficiency = *(consumption_allometric_scalar_ptr  + individuals_deme); 

			ons_constant = *(ontogenetic_niche_shift_constant_ptr  + individuals_deme); 
			ons_coefficient = *(ontogenetic_niche_shift_coefficient_ptr  + individuals_deme); 
			lw_coefficient = *(length_weight_conversion_coefficient_ptr  + individuals_deme); 
			lw_exponent = *(length_weight_conversion_exponent_ptr  + individuals_deme);
			effect_of_starvation = *(effect_of_starvation_ptr  + individuals_deme); 
			size_dependent_mortality_coefficient = *(size_dependent_mortality_coefficient_ptr  + individuals_deme);
			size_dependent_mortality_constant = *(size_dependent_mortality_constant_ptr + individuals_deme) ;

			gamma_val = *(gamma_ptr + individuals_deme);
			alpha_val = *(alpha_g_ptr + individuals_deme);

			available_resources_1 = (habitat->biotic_variables[0][individuals_deme])/(habitat->prey_array[0]->prey_maximum_abundance[individuals_deme]);
			common_satiation_resource_1 = (functional_response_scalar_resource1*available_resources_1)/(1 + handling_time_resource1 * available_resources_1); 

			available_resources_2 = (habitat->biotic_variables[1][individuals_deme])/(habitat->prey_array[1]->prey_maximum_abundance[individuals_deme]);
			common_satiation_resource_2 = (functional_response_scalar_resource2*available_resources_2)/(1 + handling_time_resource2 * available_resources_2); 

			for (int i=offset; i < offset + species->deme_sizes[j]; i++)
				{
				if (species->status[i] > 0)
					{
					satiation_resource_1 = common_satiation_resource_1 + gsl_ran_gaussian(species->gen, 0);
					satiation_resource_2 = common_satiation_resource_2 + gsl_ran_gaussian(species->gen, 0);


					xt0 = species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i];
					yt0 = species->phenotype[REVERSIBLE_MASS_PHENOTYPE][i];

					xt1 = 0;
					yt1 = 0;

					kappa = 0;
			
					maturity_threshold = 0;

					if (species->sex[i] == 0)
						maturity_threshold = F_sizes_at_maturity;
					else
						maturity_threshold = M_sizes_at_maturity;

					// Set the initial kappa value. 
					// NB: Be very careful with this: "if (x) do_this if(!x) do_that" seems to give different instructions from "if (x) do_this else do_that"
					if (xt0 >= maturity_threshold)
						kappa = 1/((1+mature_maximum_condition)*mature_maximum_condition);
					else
						kappa = 1/((1+juvenile_maximum_condition)*juvenile_maximum_condition);

					//printf("kappa: %f\n", kappa);

					// Do the dynamic energy budgeting

					maximum_consumption = 0;
					percent_adult_resource_eaten = 0;
					percent_juvenile_resource = 0;
					Eg = 0;
					prey_eaten_1 = 0;
					prey_eaten_2 = 0;

					maximum_consumption = conversion_efficiency *pow(xt0+yt0, gamma_val);

					//printf("Maximum consumption: %f\n", maximum_consumption);
					// Energy growth is therefore
					percent_adult_resource_eaten = 1/(1+exp(-(ons_constant + ons_coefficient*(lw_coefficient* pow(xt0,lw_exponent)) + gsl_ran_gaussian(species->gen, 0.0))));
					percent_juvenile_resource = 1-percent_adult_resource_eaten;

					prey_eaten_1 += satiation_resource_1*percent_juvenile_resource*maximum_consumption;
					prey_eaten_2 += satiation_resource_2*percent_adult_resource_eaten*maximum_consumption;

					//printf("Prey eaten juv. v. adult: %f %f ONS juv v. ad.: %f %f Size: %f\n", prey_eaten_1, prey_eaten_2, percent_juvenile_resource, percent_adult_resource_eaten, xt0);

					// Subtract metabolic costs
					Eg = (satiation_resource_1*percent_juvenile_resource + satiation_resource_2*percent_adult_resource_eaten)*maximum_consumption - alpha_val*(xt0+yt0);

					if (Eg >= 0)
						{
						xt1 = (yt0/xt0)*kappa*Eg + xt0;
						yt1 = (1-(yt0/xt0)*kappa)*Eg + yt0;
						}
					else
						{
						yt1 = yt0 + Eg;
						if (yt1 < 0)
							{
							yt1 = 0;
							}
						xt1=xt0;		
						}
					//printf("Juvenile Prey available: %f Energy gained: %f Original size: %f %f Final size: %f %f\n", available_resources_1, Eg, xt0, yt0, xt1, yt1);

					species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i] = xt1;
					species->phenotype[REVERSIBLE_MASS_PHENOTYPE][i] = yt1;
				
					// Update the resources eaten	
					habitat->effect_of_inds_on_biotic_variable[0][individuals_deme] += prey_eaten_1;
					habitat->effect_of_inds_on_biotic_variable[1][individuals_deme] += prey_eaten_2;


	/**********************
	Update vital rates
	**********************/
					if (Time_Step == intra_annual_time_steps - 1)
						{
						// Update fecundity for females
						if (species->sex[i] == 0)
							{
							update_female_fecundity_functor(juvenile_maximum_condition, species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i], species->phenotype[REVERSIBLE_MASS_PHENOTYPE][i], species->phenotype[EGGSIZE_PHENOTYPE][i], species->phenotype[FECUNDITY_PHENOTYPE][i]);			
							}
						}
					// Determine survivorship:
					float condition = species->phenotype[REVERSIBLE_MASS_PHENOTYPE][i] / species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i];
	
					species->phenotype[MORTALITY_PHENOTYPE][i] = (1-exp(-effect_of_starvation*condition))* (1/(1+exp(size_dependent_mortality_coefficient*(species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i] + size_dependent_mortality_constant))));

					if (species->phenotype[MORTALITY_PHENOTYPE][i] > 1)
						{
						species->phenotype[MORTALITY_PHENOTYPE][i] = 1;
						}
					if (species->phenotype[MORTALITY_PHENOTYPE][i] < 0)
						{
						species->phenotype[MORTALITY_PHENOTYPE][i] = 0;
						}

					// Leads to a maximum survivorship of about 30 years
					species->phenotype[MORTALITY_PHENOTYPE][i] = 0.9997*species->phenotype[MORTALITY_PHENOTYPE][i];
	

					rand = gsl_rng_uniform(species->gen);
				
					//printf("%f\n", species->phenotype[MORTALITY_PHENOTYPE][i] );
					//printf("Condition: %f Body size: %f Reversible mass: %f Daily Surivovrship: %f\n", condition, species->phenotype[IRREVERSIBLE_MASS_PHENOTYPE][i], species->phenotype[REVERSIBLE_MASS_PHENOTYPE][i], species->phenotype[MORTALITY_PHENOTYPE][i] );
					if (rand < species->phenotype[MORTALITY_PHENOTYPE][i])
						{
						species->status[i] = 1;
						}
					else
						{
						species->status[i] = 0;
						}
					}
				}
			}

		habitat->update();
		}
	}

void update_Fish::prepare_survivorship_constants_pointers()
	{
	effect_of_starvation_ptr = (species->demeParameters->get_vector_ptr("effect_of_starvation"));
	size_dependent_mortality_constant_ptr = (species->demeParameters->get_vector_ptr("size_dependent_mortality_constant"));
	size_dependent_mortality_coefficient_ptr = (species->demeParameters->get_vector_ptr("size_dependent_mortality_coefficient"));
	}

void update_Fish::prepare_growth_constants_pointers()
	{
	consumption_allometric_scalar_ptr = (species->demeParameters->get_vector_ptr("conversion_efficiency"));
	
	gamma_ptr = (species->demeParameters->get_vector_ptr("gamma_exponent"));
	ontogenetic_niche_shift_constant_ptr = (species->demeParameters->get_vector_ptr("ontogenetic_niche_shift_constant"));
	ontogenetic_niche_shift_coefficient_ptr = (species->demeParameters->get_vector_ptr("ontogenetic_niche_shift_coefficient"));
	length_weight_conversion_coefficient_ptr = (species->demeParameters->get_vector_ptr("length_weight_conversion_coefficient"));
	length_weight_conversion_exponent_ptr = (species->demeParameters->get_vector_ptr("length_weight_conversion_exponent"));

	handling_time_resource1_ptr = (species->demeParameters->get_vector_ptr("handling_time"));
	handling_time_resource2_ptr = (species->demeParameters->get_vector_ptr("handling_time"));

	functional_response_scalar_resource1_ptr = (species->demeParameters->get_vector_ptr("functional_response_numerator"));
	functional_response_scalar_resource2_ptr = (species->demeParameters->get_vector_ptr("functional_response_numerator"));

	/* constants governing growth rates */
	alpha_g_ptr = (species->demeParameters->get_vector_ptr("alpha_g"));

	mature_maximum_condition_ptr = (species->demeParameters->get_vector_ptr("mature_maximum_condition"));
	juvenile_maximum_condition_ptr = (species->demeParameters->get_vector_ptr("juvenile_maximum_condition"));

	M_sizes_at_maturity_ptr = (species->demeParameters->get_vector_ptr("M_sizes_at_maturity"));
	F_sizes_at_maturity_ptr = (species->demeParameters->get_vector_ptr("F_sizes_at_maturity"));

	effect_of_starvation_ptr = (species->demeParameters->get_vector_ptr("effect_of_starvation"));
	size_dependent_mortality_constant_ptr = (species->demeParameters->get_vector_ptr("size_dependent_mortality_constant"));
	size_dependent_mortality_coefficient_ptr = (species->demeParameters->get_vector_ptr("size_dependent_mortality_coefficient"));
	}

