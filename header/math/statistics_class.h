#ifndef STATISTICS_H
#define STATISTICS_H

#include <species/inds.h>
#include <species/inds_stochastic.h>
#include <util/thrust_functors.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <math/histogram.h>

#include <fstream>
#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>

class Statistics
	{
	/* 
	* A class for obtaining summary statistics about data stored in inds
	*/
	public:
		Statistics(int num_demes);
		Statistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_quantiles);
		Statistics(int num_demes, const char *output_file_quantiles);
		~Statistics();

		void calculate_mean_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		void calculate_mean_genotypes_by_deme(inds *individuals, int GENOTYPE_TO_RECORD);
		void calculate_min_max_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		
		void calculate_phenotypic_variance_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		void calculate_genotype_variance_by_deme(inds *individuals);
		void create_sample_histograms(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD);
		void calculate_quantiles(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD);

		void output_results();
		void output_histogram(int number_of_bins); // For outputting the ECDF

		void print_mean_phenotypes_by_deme();
		void print_mean_genotypes_by_deme();

	protected:
		std::ofstream summary_statistics;
		std::ofstream histogram_file;

		void genotype_phenotype_map(inds **individuals);
	
		int total_number_of_inds;
		int number_of_demes;	

		thrust::device_vector<float> mean_phenotypes;
		thrust::device_vector<float> mean_genotypes;
		thrust::device_vector<float> phenotypic_variance;
		thrust::device_vector<float> genetic_variance;
		thrust::device_vector<float> max_phenotypes;
		thrust::device_vector<float> min_phenotypes;
		thrust::device_vector<float> deme_sizes;

		thrust::device_vector<int> *histogram_by_deme;
		thrust::device_vector<float> *quantiles_by_deme;
	};


struct variance_elements_calculator
	{
	/* 
	Elements in the tuple.
	---------------------
	0: individual's random number value
	1: deme-wide mean
	2: deme size
	3: return value used to calculate variance
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		thrust::get<3>(t) = ((thrust::get<0>(t)-thrust::get<1>(t))*(thrust::get<0>(t)-thrust::get<1>(t)))/thrust::get<2>(t);
		}   
	};


#endif
