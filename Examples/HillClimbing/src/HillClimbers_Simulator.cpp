#include "HillClimbers_Simulator.h"

HillClimbers_Simulator::HillClimbers_Simulator() : Simulation()
	{
	initpop =50000*demes;
	maxpop = 50000*demes;
	
	initialize_classes();

	nspecies = 1;
	}

void HillClimbers_Simulator::initialize_classes()
	{
	array = new inds_stochastic *[nspecies];
	int species_ID = 0;

	array[0] = new HillClimbers(initpop, maxpop, seed, demes, species_ID);
	/*
	stats_hillclimber_fecundity = new Statistics(demes, "summary_statistics_fecundity.txt", 0);
	stats_hillclimber_mortality = new Statistics(demes, "summary_statistics_mortality.txt", 0);

	demographics = new DemographicStatistics(demes, "demographic_statistics.txt");

	array[0]-> exportCsv("initial_data.csv");
	*/
	}

void HillClimbers_Simulator::run()
	{
	for (int t = 0 ; t < nsteps ; t++) 
		{
		array[0]->addKids();
		array[0]->update(array);
		array[0]->removeDead();
	/*
		stats_hillclimber_fecundity->calculate_mean_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_fecundity->calculate_min_max_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_fecundity->calculate_phenotypic_variance_by_deme(array[0], array[0]->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_fecundity->output_results();

		stats_hillclimber_mortality->calculate_mean_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_mortality->calculate_min_max_phenotypes_by_deme(array[0], array[0]->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_mortality->calculate_phenotypic_variance_by_deme(array[0], array[0]->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"]);
		stats_hillclimber_mortality->output_results();

		demographics->calculate_deme_sizes(array[0]);
		demographics->record_deme_sizes();
	*/
		}
	}

HillClimbers_Simulator::~HillClimbers_Simulator()
	{
	/* cleanup */
	for (int i=0; i < nspecies; i++)
		{
		delete array[i];
		}

	delete[] array;
	/*
	delete stats_hillclimber_fecundity;
	delete stats_hillclimber_mortality;
	*/
	}
