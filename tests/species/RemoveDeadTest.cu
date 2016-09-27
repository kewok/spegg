#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <species/inds.h>

class Test_Inds_Remove_Dead : public ::testing::Test{
	protected:
		virtual void setup() {
		}
		virtual void TearDown() {
		}
};

TEST_F(Test_Inds_Remove_Dead, remove_dead) 
	{
	int initial_popsize = 10;
	int max_popsize = 20;
	int N_demes = 4;
	int species_ID = 0;
	
	inds test_inds(initial_popsize, max_popsize, N_demes, species_ID);

	// Simulate killing 
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
		test_inds.fgenotype[1][i] = i*210 + 7;

		//To try to figure out what the expected values are, based on the logic of the removeDead function:
		// std::cout << "ind phens: " << test_inds.mgenotype[0][i] << " " << test_inds.mgenotype[1][i] << "; ind deme: " << test_inds.deme[i] << "; ind status: " << test_inds.status[i] <<  std::endl;	
		}
	for (int i=initial_popsize; i  < max_popsize; i++)
		{
		test_inds.sex[i] = -242;
		test_inds.age[i] = -342;
		test_inds.id[i] = -442;
		test_inds.status[i] = -542;
		test_inds.maternal_id[i] = -642;
		test_inds.paternal_id[i] = -742;
		test_inds.deme[i] = -842;

		test_inds.phenotype[0][i] = -942;
		test_inds.phenotype[1][i] = -942;	

		test_inds.mgenotype[0][i] = -142;
		test_inds.mgenotype[1][i] = -142;

		test_inds.fgenotype[0][i] = -1042;
		test_inds.fgenotype[1][i] = -1042;
		}

	test_inds.removeDead();

	// define the expected vectors after the removeDead operation
	thrust::host_vector<int> expected_answer_status(max_popsize);
	thrust::fill(expected_answer_status.begin(), expected_answer_status.end(), 0);

	thrust::host_vector<int> expected_answer_sex(max_popsize);
	thrust::host_vector<int> expected_answer_age(max_popsize);
	thrust::host_vector<int> expected_answer_id(max_popsize);
	thrust::host_vector<int> expected_answer_maternal_id(max_popsize);
	thrust::host_vector<int> expected_answer_paternal_id(max_popsize);
	thrust::host_vector<int> expected_answer_deme(max_popsize);
	thrust::host_vector<float> expected_answer_phenotype0(max_popsize);
	thrust::host_vector<float> expected_answer_phenotype1(max_popsize);
	thrust::host_vector<float> expected_answer_mgenotype0(max_popsize);
	thrust::host_vector<float> expected_answer_mgenotype1(max_popsize);
	thrust::host_vector<float> expected_answer_fgenotype0(max_popsize);
	thrust::host_vector<float> expected_answer_fgenotype1(max_popsize);

	// Note the logic underneath remove_dead works by: (i) relabel the dead's deme as ndeme, (ii) sort the living from the dead so the living come first, (iii) among the living, sort each value by their deme.
	expected_answer_deme[0] = 1; expected_answer_deme[1] = 1; expected_answer_deme[2] = 1; 
	expected_answer_deme[3] = 3; expected_answer_deme[4] = 3; expected_answer_deme[5] = N_demes; 
	expected_answer_deme[6] = N_demes; expected_answer_deme[7] = N_demes; expected_answer_deme[8] = N_demes; expected_answer_deme[9] = N_demes;

	expected_answer_age[0] = 1; expected_answer_age[1] = 5; expected_answer_age[2] = 9; 
	expected_answer_age[3] = 3; expected_answer_age[4] = 7; expected_answer_age[5] = 0; 
	expected_answer_age[6] = 2; expected_answer_age[7] = 4; expected_answer_age[8] = 6; expected_answer_age[9] = 8;

	expected_answer_id[0] = 1; expected_answer_id[1] = 5; expected_answer_id[2] = 9; 
	expected_answer_id[3] = 3; expected_answer_id[4] = 7; expected_answer_id[5] = 0; 
	expected_answer_id[6] = 2; expected_answer_id[7] = 4; expected_answer_id[8] = 6; expected_answer_id[9] = 8;

	expected_answer_sex[0] = 1; expected_answer_sex[1] = 5; expected_answer_sex[2] = 9; 
	expected_answer_sex[3] = 3; expected_answer_sex[4] = 7; expected_answer_sex[5] = 0; 
	expected_answer_sex[6] = 2; expected_answer_sex[7] = 4; expected_answer_sex[8] = 6; expected_answer_sex[9] = 8;

	expected_answer_maternal_id[0] = 1; expected_answer_maternal_id[1] = 5; expected_answer_maternal_id[2] = 9; 
	expected_answer_maternal_id[3] = 3; expected_answer_maternal_id[4] = 7; expected_answer_maternal_id[5] = 0; 
	expected_answer_maternal_id[6] = 2; expected_answer_maternal_id[7] = 4; expected_answer_maternal_id[8] = 6; expected_answer_maternal_id[9] = 8;

	expected_answer_paternal_id[0] = 1; expected_answer_paternal_id[1] = 5; expected_answer_paternal_id[2] = 9; 
	expected_answer_paternal_id[3] = 3; expected_answer_paternal_id[4] = 7; expected_answer_paternal_id[5] = 0; 
	expected_answer_paternal_id[6] = 2; expected_answer_paternal_id[7] = 4; expected_answer_paternal_id[8] = 6; expected_answer_paternal_id[9] = 8;

	expected_answer_phenotype0[0] = 20; expected_answer_phenotype0[1] = 100; expected_answer_phenotype0[2] = 180; 
	expected_answer_phenotype0[3] = 60; expected_answer_phenotype0[4] = 140; expected_answer_phenotype0[5] = 0; 
	expected_answer_phenotype0[6] = 40; expected_answer_phenotype0[7] = 80; expected_answer_phenotype0[8] = 120; expected_answer_phenotype0[9] = 160;

	expected_answer_phenotype1[0] = 21; expected_answer_phenotype1[1] = 105; expected_answer_phenotype1[2] = 189; 
	expected_answer_phenotype1[3] = 63; expected_answer_phenotype1[4] = 147; expected_answer_phenotype1[5] = 0; 
	expected_answer_phenotype1[6] = 42; expected_answer_phenotype1[7] = 84; expected_answer_phenotype1[8] = 126; expected_answer_phenotype1[9] = 168;

	expected_answer_mgenotype0[0] = 200; expected_answer_mgenotype0[1] = 1000; expected_answer_mgenotype0[2] = 1800; 
	expected_answer_mgenotype0[3] = 600; expected_answer_mgenotype0[4] = 1400; expected_answer_mgenotype0[5] = 0; 
	expected_answer_mgenotype0[6] = 400; expected_answer_mgenotype0[7] = 800; expected_answer_mgenotype0[8] = 1200; expected_answer_mgenotype0[9] = 1600;

	expected_answer_mgenotype1[0] = 210; expected_answer_mgenotype1[1] = 1050; expected_answer_mgenotype1[2] = 1890; 
	expected_answer_mgenotype1[3] = 630; expected_answer_mgenotype1[4] = 1470; expected_answer_mgenotype1[5] = 0; 
	expected_answer_mgenotype1[6] = 420; expected_answer_mgenotype1[7] = 840; expected_answer_mgenotype1[8] = 1260; expected_answer_mgenotype1[9] = 1680;

	expected_answer_fgenotype0[0] = 207; expected_answer_fgenotype0[1] = 1007; expected_answer_fgenotype0[2] = 1807; 
	expected_answer_fgenotype0[3] = 607; expected_answer_fgenotype0[4] = 1407; expected_answer_fgenotype0[5] = 7; 
	expected_answer_fgenotype0[6] = 407; expected_answer_fgenotype0[7] = 807; expected_answer_fgenotype0[8] = 1207; expected_answer_fgenotype0[9] = 1607;

	expected_answer_fgenotype1[0] = 217; expected_answer_fgenotype1[1] = 1057; expected_answer_fgenotype1[2] = 1897; 
	expected_answer_fgenotype1[3] = 637; expected_answer_fgenotype1[4] = 1477; expected_answer_fgenotype1[5] = 7; 
	expected_answer_fgenotype1[6] = 427; expected_answer_fgenotype1[7] = 847; expected_answer_fgenotype1[8] = 1267; expected_answer_fgenotype1[9] = 1687;

	for (int i=0; i < initial_popsize; i++)
		{
		if (i < 5)
			expected_answer_status[i] = 1;
		}
	for (int i=initial_popsize; i  < max_popsize; i++)
		{
		expected_answer_sex[i] = -242;
		expected_answer_age[i] = -342;
		expected_answer_id[i] = -442;
		expected_answer_status[i] = -542;
		expected_answer_maternal_id[i] = -642;
		expected_answer_paternal_id[i] = -742;
		expected_answer_deme[i] = -842;
		expected_answer_phenotype0[i] = -942;
		expected_answer_phenotype1[i] = -942;

		expected_answer_mgenotype0[i] = -142;
		expected_answer_mgenotype1[i] = -142;

		expected_answer_fgenotype0[i] = -1042;
		expected_answer_fgenotype1[i] = -1042;
		}

	thrust::host_vector<int> host_answer_status = test_inds.status;
	thrust::host_vector<int> host_answer_deme = test_inds.deme;
	thrust::host_vector<int> host_answer_age = test_inds.age;
	thrust::host_vector<int> host_answer_id = test_inds.id;
	thrust::host_vector<int> host_answer_sex = test_inds.sex;
	thrust::host_vector<int> host_answer_paternal_id = test_inds.paternal_id;
	thrust::host_vector<int> host_answer_maternal_id = test_inds.maternal_id;
	thrust::host_vector<float> host_answer_phenotype0 = test_inds.phenotype[0];
	thrust::host_vector<float> host_answer_phenotype1 = test_inds.phenotype[1];
	thrust::host_vector<float> host_answer_mgenotype0 = test_inds.mgenotype[0];
	thrust::host_vector<float> host_answer_mgenotype1 = test_inds.mgenotype[1];
	thrust::host_vector<float> host_answer_fgenotype0 = test_inds.fgenotype[0];
	thrust::host_vector<float> host_answer_fgenotype1 = test_inds.fgenotype[1];


	EXPECT_THAT(host_answer_status, ::testing::ContainerEq(expected_answer_status));
	EXPECT_THAT(host_answer_deme, ::testing::ContainerEq(expected_answer_deme));
	EXPECT_THAT(host_answer_age, ::testing::ContainerEq(expected_answer_age));
	EXPECT_THAT(host_answer_id, ::testing::ContainerEq(expected_answer_id));
	EXPECT_THAT(host_answer_sex, ::testing::ContainerEq(expected_answer_sex));
	EXPECT_THAT(host_answer_paternal_id, ::testing::ContainerEq(expected_answer_paternal_id));
	EXPECT_THAT(host_answer_maternal_id, ::testing::ContainerEq(expected_answer_maternal_id));
	EXPECT_THAT(host_answer_phenotype0, ::testing::ContainerEq(expected_answer_phenotype0));
	EXPECT_THAT(host_answer_phenotype1, ::testing::ContainerEq(expected_answer_phenotype1));
	EXPECT_THAT(host_answer_mgenotype0, ::testing::ContainerEq(expected_answer_mgenotype0));
	EXPECT_THAT(host_answer_mgenotype1, ::testing::ContainerEq(expected_answer_mgenotype1));
	EXPECT_THAT(host_answer_fgenotype0, ::testing::ContainerEq(expected_answer_fgenotype0));
	EXPECT_THAT(host_answer_fgenotype1, ::testing::ContainerEq(expected_answer_fgenotype1));
	}

