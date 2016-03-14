#ifndef DEMOGRAPHIC_STATISTICS_H
#define DEMOGRAPHIC_STATISTICS_H

#include <math/statistics_class.h>
#include <species/inds.h>

#include <thrust/functional.h>

class DemographicStatistics : public Statistics
	{
	public:
		DemographicStatistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_histograms);
		DemographicStatistics(int num_demes, const char *output_file_summary_statistics);
		~DemographicStatistics();


		void calculate_deme_sizes(inds *species);

		void record_demographic_data();
		void record_deme_sizes();

		void calculate_age_distribution(inds *species, int number_of_bins);
		void calculate_sex_ratios(inds *species);

	protected: 
		int number_of_demes;
		thrust::device_vector<int> deme_abundances;
		thrust::device_vector<float> sex_ratios;

		void print_header();
	
	};

#endif
