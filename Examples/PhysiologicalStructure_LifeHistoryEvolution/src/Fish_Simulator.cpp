#include "Fish_Simulator.h"
#include <util/footimer2.h>

Fish_Simulator::Fish_Simulator() : Simulation()
	{
	initpop = 50000*demes;
	maxpop = 50000*demes;
	
	nspecies = 1;
	initialize_classes();
	}

void Fish_Simulator::initialize_classes()
	{
	habitat = new Fish_Habitat(seed, num_biotic_variables, num_abiotic_variables, demes, intra_step_time_steps);

	array = new inds_stochastic *[nspecies];
	int species_ID = 0;

	array[0] = new Fish(initpop, maxpop, seed + 87, demes, species_ID);

	stats_eggsize = new Statistics(demes, "summary_statistics_eggsize.txt", "useless_histograms.txt");
	stats_fecundity = new Statistics(demes, "summary_statistics_fecundity.txt", "useless_histograms2.txt");

	//array[0]-> exportCsv("initial_data.csv");

	demographics = new DemographicStatistics(demes, "demographic_statistics.txt", "age_distribution.txt");

	preyfile.open("prey_sizes.txt");
	}

void Fish_Simulator::run()
	{
	footimer2 timer, timerAll;
	timerAll.start();

	for (int t=0; t < nsteps; t++)
		{
		for (int i=0; i < nspecies; i++)
			{
			array[i]->addKids();
			
			array[i]->update(array, habitat);
			
			array[i]->removeDead();
			
			stats_eggsize->calculate_mean_phenotypes_by_deme(array[0], 5);
			stats_eggsize->calculate_min_max_phenotypes_by_deme(array[0], 5);
			stats_eggsize->calculate_phenotypic_variance_by_deme(array[0], 5);
			stats_eggsize->output_results();

			stats_fecundity->calculate_mean_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
			stats_fecundity->calculate_min_max_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
			stats_fecundity->calculate_phenotypic_variance_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
			stats_fecundity->output_results();

			demographics->calculate_deme_sizes(array[i]);
			demographics->record_deme_sizes();
			}
		habitat->update();
		
		for (int k=0; k < demes; k++)
			{
			preyfile << habitat->prey_array[0]->prey_abundance[k] << " " << habitat->prey_array[1]->prey_abundance[k] << " ";
			}
		preyfile << std::endl;
		}
	timerAll.stop();
	std::cout << "Total "; 
	timerAll.printTime();
	}

Fish_Simulator::~Fish_Simulator()
	{
	/* cleanup */
	for (int i=0; i < nspecies; i++)
		{
		delete array[i];
		}

	delete[] array;
	delete habitat;

	delete stats_fecundity;
	delete stats_eggsize;
	delete demographics;
	}

