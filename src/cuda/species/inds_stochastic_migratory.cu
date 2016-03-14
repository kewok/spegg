#include <species/inds_stochastic_migratory.h>
#include <util/file_checker.h>

inds_stochastic_migratory::inds_stochastic_migratory(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds_stochastic(size_val, maxsize_val, seed_val, ndemes, species_ID_val)
	{
/*
*
* Initialize the migration matrix. 
*
*/
	// specify the migration rate from subpopulation i to subpoulation j
	thrust::host_vector<float> migration_probability(ndemes*ndemes);

	// for convenience, suppose that the migration matrix is read in from a file Migration_Probabilities.txt
	FILE *Migration_Probabilities_File; 
	Migration_Probabilities_File = fopen("Migration_Probabilities.txt","r+");
	int available_lines = count_lines_in_file("Migration_Probabilities.txt");
	
	if (Migration_Probabilities_File == NULL)
		{
		printf("Error: specify migration probabilities!\n");
		}
	if (available_lines < ndemes*ndemes)
		{
		printf("Error: Not enough migration inputs. Simulation will not run correctly.\n");
		}
		
	for (int i=0; i < ndemes*ndemes; i++)
		{	
		float test;
		fscanf (Migration_Probabilities_File, "%f\n", &test);
		migration_probability[i] = test;
		}

	Migration_Matrix = new GSL_ProbTable *[ndemes];

	for (int i=0; i < ndemes; i++)
		{
		Migration_Matrix[i] = new GSL_ProbTable(migration_probability.begin() + i*ndemes, migration_probability.begin() + i*ndemes + ndemes);
		}
	}

inds_stochastic_migratory::~inds_stochastic_migratory()
	{
	for (int i=0; i < Num_Demes; i++)
		{
		delete Migration_Matrix[i];
		}

	delete[] Migration_Matrix;
	}
