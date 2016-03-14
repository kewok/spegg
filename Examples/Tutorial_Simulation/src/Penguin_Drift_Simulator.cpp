#include "Penguin_Drift_Simulator.h"

Penguin_Drift_Simulator::Penguin_Drift_Simulator() : Simulation()
	{
	initpop =50*demes;
	maxpop = 1000*demes;

	nspecies = 1;
	
	initialize_classes();
	}

void Penguin_Drift_Simulator::initialize_classes()
	{
	array = new inds_stochastic *[nspecies];
	int species_ID = 0;
	array[0] = new Penguins(initpop, maxpop, seed, demes, species_ID);

	stats_penguins = new Statistics(demes);

	array[0]-> exportCsv("initial_data.csv");
	}


void Penguin_Drift_Simulator::run()
	{
	for (int t = 0 ; t < nsteps ; t++) 
		{
		for (int i = 0 ; i < nspecies ; i++) 
			{
			array[i]->addKids();
			array[i]->update(array);
			array[i]->removeDead();
			}
		}

	int genotype_index_of_interest = 2;

	stats_penguins->calculate_mean_genotypes_by_deme(array[0], genotype_index_of_interest);
	std::cout << "The average genotypes at locus 1 for the two demes are:" << std::endl;
	stats_penguins->print_mean_genotypes_by_deme();

	array[0]-> exportCsv("final_data.csv");
	}


Penguin_Drift_Simulator::~Penguin_Drift_Simulator()
	{
	/* cleanup */
	for (int i=0; i < nspecies; i++)
		{
		delete array[i];
		}

	delete[] array;
	delete stats_penguins;
	}
