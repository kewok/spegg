#include <math/demographic_statistics_class.h>

DemographicStatistics::DemographicStatistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_histograms) : Statistics(num_demes, output_file_summary_statistics, output_file_histograms)
	{
	number_of_demes = num_demes;
	deme_abundances.resize(num_demes);
	sex_ratios.resize(num_demes);
	}

DemographicStatistics::DemographicStatistics(int num_demes, const char *output_file_summary_statistics) : Statistics(num_demes, output_file_summary_statistics, 0)
	{
	number_of_demes = num_demes;
	deme_abundances.resize(num_demes);
	sex_ratios.resize(num_demes);
	}

DemographicStatistics::~DemographicStatistics()
	{
	summary_statistics.close();
	if (histogram_file.is_open())
		{
		histogram_file.close();
		}
	}

void DemographicStatistics::calculate_deme_sizes(inds *species)
	{
	thrust::copy(species->deme_sizes.begin(), species->deme_sizes.begin() + number_of_demes, deme_abundances.begin());
	}

void DemographicStatistics::record_demographic_data()
	{
	for (int i=0; i < number_of_demes; i++)
		summary_statistics << deme_abundances[i] << " " << sex_ratios[i] << " " ;
	summary_statistics << std::endl;
	}


void DemographicStatistics::record_deme_sizes()
	{
	for (int i=0; i < number_of_demes; i++)
		summary_statistics << deme_abundances[i] << " ";
	summary_statistics << std::endl;
	}

void DemographicStatistics::calculate_sex_ratios(inds *species)
	{
	int ninds = species->size;
	thrust::device_vector<int> total_males(ninds);

	reduce_by_key_with_zeros(species->deme, species->sex, total_males, ninds, number_of_demes);

	calculate_deme_sizes(species);

	thrust::device_vector<float> deme_sizes(number_of_demes);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(deme_abundances.begin(), deme_sizes.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(deme_abundances.begin()+number_of_demes, deme_sizes.begin()+number_of_demes)),
		 int_to_float());

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(total_males.begin(), sex_ratios.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(total_males.begin()+number_of_demes, sex_ratios.begin()+number_of_demes)),
		 int_to_float());
	
	thrust::transform(sex_ratios.begin(), sex_ratios.begin() + number_of_demes, deme_abundances.begin(),  sex_ratios.begin(), thrust::divides<float>());
	}

void DemographicStatistics::calculate_age_distribution(inds *species, int number_of_bins)
	{
	calculate_deme_sizes(species);

	for (int i=0; i < number_of_demes; i++)
		histogram_by_deme[i].resize(number_of_bins);

	thrust::device_vector<float> ages_in_deme;

	int cumulative_deme_sizes = 0;	

	for (int i=0; i < number_of_demes; i++)
		{
		int inds_in_deme = deme_abundances[i];
		ages_in_deme.resize(inds_in_deme);

		thrust::copy(species->age.begin() + cumulative_deme_sizes, species->age.begin() + cumulative_deme_sizes + inds_in_deme, ages_in_deme.begin());
	
		if (inds_in_deme > 0)
			{
			calculate_histogram(ages_in_deme, histogram_by_deme[i], number_of_bins);
			}
		else
			{
			thrust::fill(histogram_by_deme[i].begin(), histogram_by_deme[i].begin() + number_of_bins, 0);
			}
			
		cumulative_deme_sizes += inds_in_deme; 
		}	
	}

