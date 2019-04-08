#include <math/statistics_class.h>
#include <util/gather_values_by_deme.h>

Statistics::Statistics(int num_demes)
	{
	number_of_demes = num_demes;
	
	mean_phenotypes.resize(number_of_demes);
	mean_genotypes.resize(number_of_demes);
	max_phenotypes.resize(number_of_demes);
	min_phenotypes.resize(number_of_demes);	
	phenotypic_variance.resize(number_of_demes);

	quantiles_by_deme = new thrust::device_vector<float>[number_of_demes];
	histogram_by_deme = new thrust::device_vector<int>[number_of_demes];
	}


Statistics::Statistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_histograms)
	{
	number_of_demes = num_demes;
	summary_statistics.open(output_file_summary_statistics);
	histogram_file.open(output_file_histograms);

	mean_phenotypes.resize(number_of_demes);
	mean_genotypes.resize(number_of_demes);
	max_phenotypes.resize(number_of_demes);
	min_phenotypes.resize(number_of_demes);	

	quantiles_by_deme = new thrust::device_vector<float>[number_of_demes];
	histogram_by_deme = new thrust::device_vector<int>[number_of_demes];
	}


Statistics::Statistics(int num_demes, const char *output_file_histograms)
	{
	number_of_demes = num_demes;
	histogram_file.open(output_file_histograms);

	mean_phenotypes.resize(number_of_demes);
	mean_genotypes.resize(number_of_demes);
	max_phenotypes.resize(number_of_demes);
	min_phenotypes.resize(number_of_demes);	

	quantiles_by_deme = new thrust::device_vector<float>[number_of_demes];
	histogram_by_deme = new thrust::device_vector<int>[number_of_demes];
	}


Statistics::~Statistics()
	{
	if (summary_statistics.is_open() )
		{
		summary_statistics.close();
		}

	histogram_file.close();

	delete[] histogram_by_deme;
	delete[] quantiles_by_deme;
	}

void Statistics::calculate_mean_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD)
	{
	int ninds = individuals->size;

	reduce_by_key_with_zeros(individuals->deme, individuals->phenotype[PHENOTYPE_TO_RECORD], mean_phenotypes, ninds, number_of_demes);

	deme_sizes.resize(number_of_demes);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals->deme_sizes.begin(), deme_sizes.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(individuals->deme_sizes.begin()+number_of_demes, deme_sizes.begin()+number_of_demes)),
		 int_to_float());
	
	thrust::transform(mean_phenotypes.begin(), mean_phenotypes.begin() + number_of_demes, deme_sizes.begin(),  mean_phenotypes.begin(), thrust::divides<float>());
	}

void Statistics::calculate_mean_genotypes_by_deme(inds *individuals, int GENOTYPE_TO_RECORD)
	{
	int ninds = individuals->size;
	thrust::device_vector<float> summed_genotype(ninds);
	
	thrust::transform(individuals->fgenotype[GENOTYPE_TO_RECORD].begin(), individuals->fgenotype[GENOTYPE_TO_RECORD].begin() + ninds, individuals->mgenotype[GENOTYPE_TO_RECORD].begin(), summed_genotype.begin(), thrust::plus<float>()); 
	thrust::device_vector<float> twos(ninds);
	thrust::fill(twos.begin(), twos.end(), 2.0);
	thrust::transform(summed_genotype.begin(), summed_genotype.begin() + ninds, twos.begin(), summed_genotype.begin(),  thrust::divides<float>());


	reduce_by_key_with_zeros(individuals->deme, summed_genotype, mean_genotypes, ninds, number_of_demes);

	deme_sizes.resize(number_of_demes);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals->deme_sizes.begin(), deme_sizes.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(individuals->deme_sizes.begin()+number_of_demes, deme_sizes.begin()+number_of_demes)),
		 int_to_float());
	
	thrust::transform(mean_genotypes.begin(), mean_genotypes.begin() + number_of_demes, deme_sizes.begin(),  mean_genotypes.begin(), thrust::divides<float>());
	}

void Statistics::calculate_min_max_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD)
	{
	int cumulative_deme_sizes = 0;
	for (int i=0; i < number_of_demes; i++)
		{
		if (individuals->deme_sizes[i] > 0)
			{
			max_phenotypes[i] = *thrust::max_element(individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes, individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes + individuals->deme_sizes[i]);
			min_phenotypes[i] = *thrust::min_element(individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes, individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes + individuals->deme_sizes[i]);
			}
		else
			{
			max_phenotypes[i] = 0;
			min_phenotypes[i] = 0;
			}
		cumulative_deme_sizes += individuals->deme_sizes[i];
		}
	}

void Statistics::calculate_phenotypic_variance_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD)
	{
	// Note: the current code assumes calculate_mean_phenotypes_by_deme() has always been run just prior to this code. I'm not sure that is a safe assumption, and we might need to modify that to perform the calculations of means inside calculate_phenotypic_variance_by_deme.

	int ninds = individuals->size;

	thrust::device_vector<float> squared_phenotype_values(ninds);
	thrust::device_vector<float> individualized_mean_phenotypes(ninds);
	thrust::device_vector<float> individualized_deme_sizes(ninds);
	thrust::device_vector<int> individual_indices(ninds);
	thrust::sequence(individual_indices.begin(), individual_indices.begin() + individuals->size, 0);

	gather_values_by_deme(individual_indices, individuals->deme, mean_phenotypes, individualized_mean_phenotypes);
	gather_values_by_deme(individual_indices, individuals->deme, deme_sizes, individualized_deme_sizes);

	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals->phenotype[PHENOTYPE_TO_RECORD].begin(), individualized_mean_phenotypes.begin(), individualized_deme_sizes.begin(), squared_phenotype_values.begin())),
		 thrust::make_zip_iterator(thrust::make_tuple(individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + ninds, individualized_deme_sizes.begin()+ninds, individualized_mean_phenotypes.begin()+ninds, squared_phenotype_values.begin() + ninds)),
		 variance_elements_calculator());

	reduce_by_key_with_zeros(individuals->deme, squared_phenotype_values, phenotypic_variance, ninds, number_of_demes);
	}

void Statistics::output_results()
	{
	// If you aren't tracking the variances:
	if (phenotypic_variance.size() == 0)
		{
		phenotypic_variance.resize(number_of_demes);
		}

	for (int i=0; i < number_of_demes; i++)
		{
		summary_statistics << min_phenotypes[i] << " " << mean_phenotypes[i] << " " << max_phenotypes[i] << " " << phenotypic_variance[i] << " ";
		}
	summary_statistics << std::endl;
	}

void Statistics::output_histogram(int number_of_bins)
	{
/*	for (int i=0; i < number_of_demes; i++)
		{
		for (int j=0; j < number_of_bins; j++)
			summary_statistics << quantiles_by_deme[i][j] << " ";
		}
*/	
	for (int i=0; i < number_of_demes; i++)
		{
		for (int j=0; j < number_of_bins; j++)
			{
			histogram_file << histogram_by_deme[i][j] << " ";
			}
		}
	histogram_file << std::endl;
	}


void Statistics::calculate_quantiles(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD)
	{
	for (int i=0; i < number_of_demes; i++)
		quantiles_by_deme[i].resize(number_of_bins);

	thrust::device_vector<float> phenotypes_in_deme;

	int cumulative_deme_sizes = 0;	

	for (int i=0; i < number_of_demes; i++)
		{
		int inds_in_deme = individuals->deme_sizes[i];
		phenotypes_in_deme.resize(inds_in_deme);

		thrust::device_vector<int> draw_at(number_of_bins);

		int interval = (int) inds_in_deme/number_of_bins;

		thrust::sequence(draw_at.begin(), draw_at.end(), 0, interval);

		if (i==0)
			{
			thrust::copy(individuals->phenotype[PHENOTYPE_TO_RECORD].begin(), individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + inds_in_deme, phenotypes_in_deme.begin());
			}
		else
			{
			thrust::copy(individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes, individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes + inds_in_deme, phenotypes_in_deme.begin());
			}

		thrust::sort(phenotypes_in_deme.begin(), phenotypes_in_deme.begin() + inds_in_deme);
		
		thrust::gather(draw_at.begin(), draw_at.end(), phenotypes_in_deme.begin(), quantiles_by_deme[i].begin());


		cumulative_deme_sizes += inds_in_deme; 
		}	
	}

void Statistics::create_sample_histograms(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD)
	{
	for (int i=0; i < number_of_demes; i++)
		histogram_by_deme[i].resize(number_of_bins);

	thrust::device_vector<float> phenotypes_in_deme;

	int cumulative_deme_sizes = 0;	

	for (int i=0; i < number_of_demes; i++)
		{
		int inds_in_deme = individuals->deme_sizes[i];
		phenotypes_in_deme.resize(inds_in_deme);

		thrust::copy(individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes, individuals->phenotype[PHENOTYPE_TO_RECORD].begin() + cumulative_deme_sizes + inds_in_deme, phenotypes_in_deme.begin());
	
		if (inds_in_deme > 0)
			{
			calculate_histogram(phenotypes_in_deme, histogram_by_deme[i], number_of_bins);
			}
		else
			{
			thrust::fill(histogram_by_deme[i].begin(), histogram_by_deme[i].begin() + number_of_bins, 0);
			}
			
		cumulative_deme_sizes += inds_in_deme; 
		}	
	}

void Statistics::print_mean_phenotypes_by_deme()
	{
	for (int i=0; i < mean_phenotypes.size(); i++)
		std::cout << mean_phenotypes[i] << " ";
	std::cout << std::endl;
	}

void Statistics::print_mean_genotypes_by_deme()
	{
	for (int i=0; i < mean_genotypes.size(); i++)
		std::cout << mean_genotypes[i] << " ";
	std::cout << std::endl;
	}
