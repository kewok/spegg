#include "coevolutionSimulator.h"

Coevolution_Simulator::Coevolution_Simulator() : Simulation()
	{
	initpop = 25000*demes;
	maxpop = 100000*demes;
	
	nspecies = 2;
	initialize_classes();
	}

void Coevolution_Simulator::initialize_classes()
	{
	array = new inds_stochastic *[nspecies];
	for (int i=0; i < nspecies; i++)
		{
		array[i] = initialize_species(initpop, maxpop, seed, demes, i);
		}
	
	stats_prey_competitive_ability = new Statistics(demes, "summary_statistics_resource_competitive_ability.txt", "histograms_competitive_ability_frequency.txt");
	stats_predator_attack_ability = new Statistics(demes, "summary_statistics_consumer_attack_ability.txt", "histograms_consumer_attack_frequency.txt");
	stats_prey_defense_ability = new Statistics(demes, "summary_statistics_resource_defense.txt", "histograms_resource_defense_frequency.txt");;
	demographics_predator = new DemographicStatistics(demes, "demographic_statistics_consumer.txt");	
	demographics_prey = new DemographicStatistics(demes, "demographic_statistics_resource.txt");
	
	}

void Coevolution_Simulator::run()
	{
	for (int t=0; t < nsteps; t++)
		{
		for (int i=0; i < nspecies; i++)
			{
			array[i]->addKids();
			array[i]->update(array);

			// Since array[i] is a pointer to an inds_stochastic object, you have to convert it to a pointer for the derived inds_stochastic_migratory class via static_cast<derived *> in order for the compilation to work. 
			static_cast<inds_stochastic_migratory *> (array[i])->migrate();

			array[i]->removeDead();

			if (i==0)
					{
					stats_predator_attack_ability->calculate_mean_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["CONSUMER_ATTACK_PHENOTYPE_INDEX"]);
					stats_predator_attack_ability->calculate_min_max_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["CONSUMER_ATTACK_PHENOTYPE_INDEX"]);
					stats_predator_attack_ability->calculate_phenotypic_variance_by_deme(array[i], array[i]->demeParameters->species_specific_values["CONSUMER_ATTACK_PHENOTYPE_INDEX"]);
					stats_predator_attack_ability->output_results();
					demographics_predator->calculate_deme_sizes(array[i]);
					demographics_predator->record_deme_sizes();
					}
			if (i==1)
					{
					stats_prey_defense_ability->calculate_mean_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["RESOURCE_DEFENSE_PHENOTYPE_INDEX"]);
					stats_prey_defense_ability->calculate_min_max_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["RESOURCE_DEFENSE_PHENOTYPE_INDEX"]);
					stats_prey_defense_ability->calculate_phenotypic_variance_by_deme(array[i], array[i]->demeParameters->species_specific_values["RESOURCE_DEFENSE_PHENOTYPE_INDEX"]);
					stats_prey_defense_ability->output_results();

					stats_prey_competitive_ability->calculate_mean_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["RESOURCE_COMPETITION_PHENOTYPE_INDEX"]);
					stats_prey_competitive_ability->calculate_min_max_phenotypes_by_deme(array[i], array[i]->demeParameters->species_specific_values["RESOURCE_COMPETITION_PHENOTYPE_INDEX"]);
					stats_prey_competitive_ability->output_results();

					demographics_prey->calculate_deme_sizes(array[i]);
					demographics_prey->record_deme_sizes();
					}
			}
		}
	}

Coevolution_Simulator::~Coevolution_Simulator()
	{
	for (int i=0; i < nspecies; i++)
		{
		delete array[i];
		}
	delete[] array;
	}

coevolvingSpecie *Coevolution_Simulator::initialize_species(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val)
	{
	// Arbitrary designate species 0 to be the consumer and species 1 to be the resource
	if (species_ID_val == 0)
		return new myConsumer(size_val, maxsize_val, seed_val, ndemes, species_ID_val);

	if (species_ID_val == 1)
		return new myResource(size_val, maxsize_val, seed_val, ndemes, species_ID_val);
	}
