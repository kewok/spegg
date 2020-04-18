#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <species/inds.h>

class Test_Inds_Initialize_CSV : public ::testing::Test{
	protected:
		virtual void setup() {
		}
		virtual void TearDown() {
		}
};

TEST_F(Test_Inds_Initialize_CSV, initialize_csv) 
	{
	int initial_popsize = 10;
	int max_popsize = 20;
	int N_demes = 4;
	int species_ID = 0;
	
	inds test_inds(initial_popsize, max_popsize, N_demes, species_ID);

	// Do a manual initialization
	for (int i=0; i < initial_popsize; i++)
		{
		test_inds.id[i] = i;
		test_inds.age[i] = i;
		test_inds.sex[i] = i;
		test_inds.maternal_id[i] = i;
		test_inds.paternal_id[i] = i;
		test_inds.deme[i] = i % 4;
		test_inds.status[i] = i % 2;

		test_inds.phenotype[0][i] = i*20;
		test_inds.phenotype[1][i] = i*21;

		test_inds.mgenotype[0][i] = i*200;
		test_inds.mgenotype[1][i] = i*210;

		test_inds.fgenotype[0][i] = i*200 + 7;
		test_inds.fgenotype[1][i] = i*210 + 7;;	
		}

	// define the expected vectors after the modification and re-initialization operation
	thrust::host_vector<int> expected_answer_id(max_popsize);
	thrust::copy(test_inds.id.begin(), test_inds.id.begin() + max_popsize,  expected_answer_id.begin());
	thrust::host_vector<float> expected_answer_phenotype1(max_popsize);
	thrust::copy(test_inds.phenotype[1].begin(), test_inds.phenotype[1].begin() + max_popsize,  expected_answer_phenotype1.begin());
	thrust::host_vector<float> expected_answer_mgenotype0(max_popsize);
	thrust::copy(test_inds.mgenotype[0].begin(), test_inds.mgenotype[0].begin() + max_popsize,  expected_answer_mgenotype0.begin());
	thrust::host_vector<float> expected_answer_fgenotype1(max_popsize);
	thrust::copy(test_inds.fgenotype[1].begin(), test_inds.fgenotype[1].begin() + max_popsize,  expected_answer_fgenotype1.begin());

	test_inds.exportCsv("test_file.csv");

	// Make some arbitrary modifications
	for (int i=0; i < initial_popsize; i++)
		{
		test_inds.id[i] = 0;
		test_inds.phenotype[1][i] = i*748;
		test_inds.mgenotype[0][i] = i*52;
		test_inds.fgenotype[1][i] = i*53;
		}

	thrust::host_vector<int> host_answer_id = test_inds.id;
	thrust::host_vector<float> host_answer_phenotype1 = test_inds.phenotype[1];
	thrust::host_vector<float> host_answer_mgenotype0 = test_inds.mgenotype[0];
	thrust::host_vector<float> host_answer_fgenotype1 = test_inds.fgenotype[1];

	// First, confirm the modification changed the values:
	std::cout << " Check that modification occured" << std::endl;
	EXPECT_THAT(host_answer_id, ::testing::Ne(expected_answer_id));
	EXPECT_THAT(host_answer_phenotype1, ::testing::Ne(expected_answer_phenotype1));
	EXPECT_THAT(host_answer_mgenotype0, ::testing::Ne(expected_answer_mgenotype0));
	EXPECT_THAT(host_answer_fgenotype1, ::testing::Ne(expected_answer_fgenotype1));

	test_inds.initialize_from_CSV("test_file.csv");

	// Check that initialization from CSV file works:
	thrust::host_vector<int> host_answer_id_revised = test_inds.id;
	thrust::host_vector<float> host_answer_phenotype1_revised = test_inds.phenotype[1];
	thrust::host_vector<float> host_answer_mgenotype0_revised = test_inds.mgenotype[0];
	thrust::host_vector<float> host_answer_fgenotype1_revised = test_inds.fgenotype[1];

	EXPECT_THAT(host_answer_id_revised, ::testing::ContainerEq(expected_answer_id));
	EXPECT_THAT(host_answer_phenotype1_revised, ::testing::ContainerEq(expected_answer_phenotype1));
	EXPECT_THAT(host_answer_mgenotype0_revised, ::testing::ContainerEq(expected_answer_mgenotype0));
	EXPECT_THAT(host_answer_fgenotype1_revised, ::testing::ContainerEq(expected_answer_fgenotype1));


	remove("test_file.csv");
	}

