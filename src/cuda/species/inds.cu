#include <species/inds.h>
#include <math/thrust_probabilities.h>

#include <util/rapidcsv/src/rapidcsv.h>

#include <curand.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <thrust/count.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/scatter.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>


inds::inds(int size_val, int maxsize_val, int num_demes, int species_ID_val) : size(size_val), maxsize(maxsize_val), Num_Demes(num_demes), nextid(size_val), species_ID(species_ID_val)
	{
	/*
	* A very limited initialization method which creates the data structures and performs a basic sanity check to make sure that the maximum and starting numbers of individuals are biologically meaningful.
	*/
	//Sanity check.
	if (size_val < 0 || maxsize_val < 0) {
		std::cerr << "Population size must be non-negative!" << std::endl;
		exit(1);
	}
	if (size_val > maxsize_val) {
		std::cerr << "Maximum size must be greater or equal to initial population size!" << std::endl;
		exit(1);
	}
	deme_sizes.resize(Num_Demes);
	max_deme_sizes.resize(Num_Demes);

	demeParameters = new DemeSettings("deme_config.txt", species_ID);

	if (demeParameters->check_number_of_demes() < Num_Demes)
		{
		std::cout << "The number of demes specified in the simulation_conf file cannot exceed the number specified in your deme_config.txt file; Please fix this before proceeding. " << std::endl;
		exit(1);
		}

	nloci = (int) demeParameters->GeneticArchitecture->Number_of_Loci;
	nphen = (int) demeParameters->GeneticArchitecture->Number_of_Phenotypes;

	initialize_individuals(nloci, nphen);
	//Set maxsize.
	setMaxSize(maxsize);

	demeCalculations();
	//Fill in ID, STATUS, and DEME.
	thrust::sequence(id.begin(), id.begin() + size);
	thrust::fill(status.begin(), status.begin() + size, 1);
	}

void inds::initialize_individuals(int nloci, int nphen)
	{
	//Allocate gene and phen data vectors.
	fgenotype = new thrust::device_vector<float>[nloci];
	mgenotype = new thrust::device_vector<float>[nloci];
	phenotype = new thrust::device_vector<float>[nphen];
	}

void inds::initialize_from_CSV(const char *filename)
	{
	rapidcsv::Document csv(filename);

	// Check that MaxSize is satisfied and exit if it isn't. 
	size = csv.GetRowCount();
	if (size > maxsize)
		{
		std::cerr << "Maximum size must be greater or equal to initial population size!" << std::endl;
		exit(1);
		}

	std::vector<int> id_in = csv.GetColumn<int>("Id");
	thrust::copy(id_in.begin(), id_in.end(), id.begin());

	std::vector<int> status_in = csv.GetColumn<int>("Status");
	thrust::copy(status_in.begin(), status_in.end(), status.begin());

	std::vector<int> sex_in = csv.GetColumn<int>("Sex");
	thrust::copy(sex_in.begin(), sex_in.end(), sex.begin());

	std::vector<int> age_in = csv.GetColumn<int>("Age");
	thrust::copy(age_in.begin(), age_in.end(), age.begin());

	std::vector<int> deme_in = csv.GetColumn<int>("Deme");
	thrust::copy(deme_in.begin(), deme_in.end(), deme.begin());
	
	std::string mgen_in = "mgene";
	std::string fgen_in = "fgene";

	for (int i=0; i < nloci; i++)
		{
		std::stringstream tempm;
		tempm << mgen_in << i ;
		std::stringstream tempf;
		tempf << fgen_in << i ;
	
		std::vector<float> mgenotype_in = csv.GetColumn<float>(tempm.str());
		std::vector<float> fgenotype_in = csv.GetColumn<float>(tempf.str());
		thrust::copy(mgenotype_in.begin(), mgenotype_in.end(), mgenotype[i].begin());
		thrust::copy(fgenotype_in.begin(), fgenotype_in.end(), fgenotype[i].begin());
		}

	std::string phen_in = "phen";
	for (int i=0; i < nphen; i++)
		{
		std::stringstream tempPh;
		tempPh << phen_in << i;
	
		std::vector<float> phenotype_in = csv.GetColumn<float>(tempPh.str());
		thrust::copy(phenotype_in.begin(), phenotype_in.end(), phenotype[i].begin());
		}

	// Optionally import parental IDs
	std::vector<std::string> columnNames = csv.GetColumnNames();
	bool maternal_ID_exists = (std::find(columnNames.begin(), columnNames.end(), "maternal_ID") != columnNames.end());
	if (maternal_ID_exists)
		{
		std::vector<int> maternal_id_in = csv.GetColumn<int>("maternal_ID");
		thrust::copy(maternal_id_in.begin(), maternal_id_in.end(), maternal_id.begin());
		}

	bool paternal_ID_exists = (std::find(columnNames.begin(), columnNames.end(), "paternal_ID") != columnNames.end());
	if (paternal_ID_exists)
		{
		std::vector<int> paternal_id_in = csv.GetColumn<int>("paternal_ID");
		thrust::copy(paternal_id_in.begin(), paternal_id_in.end(), paternal_id.begin());
		}
	}

inds::~inds()
	{
	/*
	* The destructor to clear GPU RAM; . Note only the thrust vectors of arrays are destroyed.
	*/
	delete[] fgenotype;
	delete[] mgenotype;
	delete[] phenotype;
	}

void inds::exportCsv()
	{
	std::ostringstream stream;
	stream << "species_" << species_ID << ".csv";
	std::string result = stream.str();
	const char *mySpeciesTitle = result.c_str();
	exportCsv(mySpeciesTitle);
	}

void inds::exportCsv(int timestep)
	{
	std::ostringstream stream;
	stream << "species_" << species_ID << "_" << timestep << ".csv";
	std::string result = stream.str();
	const char *mySpeciesTitle = result.c_str();
	exportCsv(mySpeciesTitle);
	}

void inds::exportCsv(const char *filename)
	{
	/*
	* Migrate all species-level data stored on the GPU ram to a CSV file on the main drive at any arbirary point in time.
	*/
	std::ofstream file;
	file.open(filename);
	
	//Output size
	//file << "Size,NLoci,NPhen,Seed" << std::endl;
	//file << size << "," << nloci << "," << nphen << "," << seed << std::endl << std::endl;
	
	//Output headers
	file << "Index,Id,Status,Sex,Age,Deme,";
	for (int i = 0 ; i < nloci ; i++) {
		file << "fgene" << i << ",";
	}
	for (int i = 0 ; i < nloci ; i++) {
		file << "mgene" << i << ",";
	}

	for (int i = 0 ; i < nphen ; i++) {
		if (i < nphen - 1)
			file << "phen" << i << ",";
		else
			file << "phen" << i << std::endl;
	}
	
	//Output n individuals
	for (int i = 0 ; i < size ; i++) {
		file << i << "," << id[i] << "," << status[i] << "," << sex[i] << "," << age[i] << "," << deme[i] << ",";
		for (int j = 0 ; j < nloci ; j++) {
			file << fgenotype[j][i] << ",";
		}
		for (int j = 0 ; j < nloci ; j++) {
			file << mgenotype[j][i] << ",";
		}
	for (int j = 0 ; j < nphen ; j++) {
		if (j < nphen - 1)
			file << phenotype[j][i] << ",";
		else
			file << phenotype[j][i] << std::endl;
		}
	}
	
	file.close();
	}

void inds::exportCsv(const char *filename, int timestep)
	{
	/*
	* An overloaded version of the export CSV file which prints parental ids of the individuals, as well as the time step yr/in which their data are recorded to the CSV file.
	*/
	std::ofstream file;
	file.open(filename,std::ios_base::app);
	
	//Output size
	//file << "Size,NLoci,NPhen,Seed" << std::endl;
	//file << size << "," << nloci << "," << nphen << "," << seed << std::endl << std::endl;
	
	//Output headers
	if ((timestep==0))
		{
		file << "Time_step, Index,Id,Status,Sex,Age,Deme,maternal_ID,paternal_ID,";
		for (int i = 0 ; i < nloci ; i++) {
			file << "fgene" << i << ",";
		}
		for (int i = 0 ; i < nloci ; i++) {
			file << "mgene" << i << ",";
		}
	
		for (int i = 0 ; i < nphen ; i++) {
			if (i < nphen - 1)
				file << "phen" << i << ",";
			else
				file << "phen" << i << std::endl;
		}
		}
	
	//Output n individuals
	for (int i = 0 ; i < size ; i++) {
		file << timestep << "," << i << "," << id[i] << "," << status[i] << "," << sex[i] << "," << age[i] << "," << deme[i] << "," << maternal_id[i] << "," << paternal_id[i] << ",";
		for (int j = 0 ; j < nloci ; j++) {
			file << fgenotype[j][i] << ",";
		}
		for (int j = 0 ; j < nloci ; j++) {
			file << mgenotype[j][i] << ",";
		}
	
		for (int j = 0 ; j < nphen ; j++) {
			if (j < nphen - 1)
				file << phenotype[j][i] << ",";
			else
				file << phenotype[j][i] << std::endl;
		}
	}
	
	file.close();
	}

void inds::exportCsv(const char *filename, int timestep1, int timestep2)
	{
	/*
	*Output file stream for each time-step timestep within yr which is incremented in main.cpp
	*/
	std::ofstream file;
	file.open(filename,std::ios_base::app);
	
	//Output size
	//file << "Size,NLoci,NPhen,Seed" << std::endl;
	//file << size << "," << nloci << "," << nphen << "," << seed << std::endl << std::endl;
	
	//Output headers
	if ((timestep1==0) && (timestep2 == 0))
		{
		file << "Time_Step1,Time_Step2,Index,Id,Status,Sex,Age,Deme,maternal_ID,paternal_ID";
	/*	for (int i = 0 ; i < nloci ; i++) {
			file << "fgene" << i << ",";
		}
		for (int i = 0 ; i < nloci ; i++) {
			file << "mgene" << i << ",";
		}
	*/
		for (int i = 0 ; i < nphen ; i++) {
			if (i < nphen - 1)
				file << "phen" << i << ",";
			else
				file << "phen" << i << std::endl;
			}
		}
	
	//Output n individuals
	for (int i = 0 ; i < size ; i++) {
		file << timestep1 << "," << timestep2 << "," << i << "," << id[i] << "," << status[i] << "," << sex[i] << "," << age[i] << "," << deme[i] << "," << maternal_id[i] << "," << paternal_id[i] << ",";
	/*	for (int j = 0 ; j < nloci ; j++) {
			file << fgenotype[j][i] << ",";
		}
		for (int j = 0 ; j < nloci ; j++) {
			file << mgenotype[j][i] << ",";
		}
	*/
	for (int j = 0 ; j < nphen ; j++) {
		if (j < nphen - 1)
			file << phenotype[j][i] << ",";
		else
			file << phenotype[j][i] << std::endl;
		}
	}
	
	file.close();
	}

void inds::removeDead()
	{
	/*
	* A function to reorganize your data points so that inds (0, 1, 2, ..., size) consist only of individuals whose vital status = 1.
	* This method rearranges the data structures of inds so that all the data points represented in each data structure from 0 to number_of_individuals represent individuals that are alive. The data points past number_of_individuals are garbage. All operations on inds or its derived classes should, therefore, operate only on number_of_individuals data points. As a rule, number_of_individuals <= max_number_of_individuals; when the two are equal, the data from dead individuals should simply be overwritten until number_of_individuals < max_number_of_individuals again.

	* Ideally this would involve stream compaction rather than relegating the dead individuals to occupy empty spaces, but some preliminary experiments suggested stream compaction entails substantial performance costs in thrust compared to the approach using sorting and gathering. This solution is subject to change; should performance improvements in Thrust or CUDA allow it, we will return to stream compaction routines here instead.
	*
	*/
	//Keys to flag those to be removed.
	thrust::device_vector<int> keys(size);
	
	//Setup keys. 1 = to be removed. 0 = not removed.
	thrust::transform(status.begin(), status.begin() + size, keys.begin(), thrust::logical_not<int>());

	reassign_dead_deme_functor give_deads_deme(Num_Demes);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(status.begin(), deme.begin())),
        		 thrust::make_zip_iterator(thrust::make_tuple(status.begin() + size, deme.begin() + size)),
			give_deads_deme);

	//Count number of alive individuals.
	int num_alive = thrust::count(status.begin(), status.begin() + size, 1);
	
	// Perform the sort based on deme with max_demes + 1 at the bottom
	// Then perform gather
	
	thrust::device_vector<int> map(size);
	thrust::copy(deme.begin(), deme.begin() + size, keys.begin());
	thrust::sequence(map.begin(), map.begin() + size);
	thrust::stable_sort_by_key(keys.begin(), keys.begin() + size, map.begin(), thrust::less<int>());

	thrust::device_vector<int> new_vals(size);
	thrust::device_vector<float> new_vals_float(size);

	thrust::gather(map.begin(), map.begin() + size, id.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, id.begin());

	thrust::gather(map.begin(), map.begin() + size, age.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, age.begin());


	thrust::gather(map.begin(), map.begin() + size, sex.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, sex.begin());
	
	thrust::gather(map.begin(), map.begin() + size, status.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, status.begin());

	thrust::gather(map.begin(), map.begin() + size, maternal_id.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, maternal_id.begin());

	thrust::gather(map.begin(), map.begin() + size, paternal_id.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, paternal_id.begin());

	thrust::gather(map.begin(), map.begin() + size, deme.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + size, deme.begin());

	for (int i = 0 ; i < nloci ; i++) {
		thrust::gather(map.begin(), map.begin() + size, fgenotype[i].begin(), new_vals_float.begin());
		thrust::copy(new_vals_float.begin(), new_vals_float.begin() + size, fgenotype[i].begin());
		thrust::gather(map.begin(), map.begin() + size, mgenotype[i].begin(), new_vals_float.begin());
		thrust::copy(new_vals_float.begin(), new_vals_float.begin() + size, mgenotype[i].begin());
		}

	for (int i = 0 ; i < nphen ; i++) {
		thrust::gather(map.begin(), map.begin() + size, phenotype[i].begin(), new_vals_float.begin());
		thrust::copy(new_vals_float.begin(), new_vals_float.begin() + size, phenotype[i].begin());
		}

	size = num_alive;
	// Recalculate the number of individuals per deme
	demeCalculations();
	}

void inds::setMaxSize(int n)
	{
	/*
	* Allocate memory RAM on the GPU device to be used in your inds class during sPEGG simulation. 
	*/
	if (n < 0) {
		std::cerr << "setMaxSize: size must be non-negative!" << std::endl;
		return;
	}
	
	id.resize(n);
	status.resize(n);
	sex.resize(n);
	age.resize(n);
	deme.resize(n);

	// Initialize vectors storing parental arrays:
	maternal_id.resize(n);
	thrust::fill(maternal_id.begin(), maternal_id.begin() + n, -1);

	paternal_id.resize(n);
	thrust::fill(paternal_id.begin(), paternal_id.begin() + n, -1);
	
	for (int i = 0 ; i < nloci ; i++) {
		fgenotype[i].resize(n);
		mgenotype[i].resize(n);
	}
	for (int i = 0 ; i < nphen ; i++) {
		phenotype[i].resize(n);
	}
	
	maxsize = n;
	}

void inds::sortByDeme()
	{
	/* 
	* Reorganize the individual-level data points as needed, according to the deme to which the individual belongs.  As in \link removeDead() 
the data points past number_of_individuals are garbage. Note that the vital status for all individuals in sort_by_deme is assumed to be 1.
	*
	*/
	//Count number of alive individuals.
	int num_alive = thrust::count(status.begin(), status.begin() + size, 1);

	//Keys
	thrust::device_vector<int> keys(num_alive);
	
	/*
		For each vector,
		set keys to the deme vector.
		Then use the keys to sort.
	*/

	// Relocate the values according to the map generated by sort
	thrust::device_vector<int> map(num_alive);
	thrust::copy(deme.begin(), deme.begin() + num_alive, keys.begin());
	thrust::sequence(map.begin(), map.begin() + num_alive);
	thrust::stable_sort_by_key(keys.begin(), keys.begin() + num_alive, map.begin());

	// Map now directs where individuals should go

	thrust::device_vector<int> new_vals(size);
	thrust::device_vector<float> new_vals_float(size);
	
	thrust::gather(map.begin(), map.begin() + num_alive, id.begin(), new_vals.begin());

	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, id.begin());
	thrust::gather(map.begin(), map.begin() + num_alive, age.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, age.begin());
	

	thrust::gather(map.begin(), map.begin() + num_alive, sex.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, sex.begin());
	
	thrust::gather(map.begin(), map.begin() + num_alive, status.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, status.begin());

	thrust::gather(map.begin(), map.begin() + num_alive, maternal_id.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, maternal_id.begin());

	thrust::gather(map.begin(), map.begin() + num_alive, paternal_id.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, paternal_id.begin());

	thrust::gather(map.begin(), map.begin() + num_alive, deme.begin(), new_vals.begin());
	thrust::copy(new_vals.begin(), new_vals.begin() + num_alive, deme.begin());

	for (int i = 0 ; i < nloci ; i++) 
		{
		thrust::gather(map.begin(), map.begin() + num_alive, fgenotype[i].begin(), new_vals_float.begin());
		thrust::copy(new_vals_float.begin(), new_vals_float.begin() + num_alive, fgenotype[i].begin());
		thrust::gather(map.begin(), map.begin() + num_alive, mgenotype[i].begin(), new_vals_float.begin());
		thrust::copy(new_vals_float.begin(), new_vals_float.begin() + num_alive, mgenotype[i].begin());
		}

	for (int i = 0 ; i < nphen ; i++) {
			thrust::gather(map.begin(), map.begin() + num_alive, phenotype[i].begin(), new_vals_float.begin());
			thrust::copy(new_vals_float.begin(), new_vals_float.begin() + num_alive, phenotype[i].begin());
		}
	}

void inds::demeCalculations()
	{
	/*
	* Determine how many individuals are in each deme. The results are stored in vector deme_sizes, which is accessed via get_deme_sizes(thrust::device_vector<int> &deme_sizevals) where the vector deme_sizevals has already been declared and instantiated.
	*/
	thrust::counting_iterator<int> search_begin(0);
	thrust::device_vector<int> temp_deme_sizes;

	temp_deme_sizes.resize(Num_Demes);

	thrust::upper_bound(deme.begin(), deme.begin() + size,
                      search_begin, search_begin + Num_Demes,
                      temp_deme_sizes.begin());

	thrust::adjacent_difference(temp_deme_sizes.begin(), temp_deme_sizes.end(),
                              deme_sizes.begin());
	}
