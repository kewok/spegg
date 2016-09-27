#include <util/gather_values_by_deme.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class GatherValuesByDemeIntSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(GatherValuesByDemeIntSuccess, gather_values_by_deme) {
	thrust::device_vector<int> indices(4);
	thrust::device_vector<int> demes(10);
	thrust::device_vector<int> deme_specific_values(3);
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(4);

	indices[0] = 0; indices[1] = 5; indices[2] = 8; indices[3] = 2;

	demes[0] = 0; demes[1] = 0; demes[2] = 1; demes[3] = 2; demes[4] = 0; demes[5] = 2; demes[6] = 1; demes[7]=2; demes[8]=1;  demes[9] = 0;
	
	deme_specific_values[0] = 100; deme_specific_values[1] = 1000; deme_specific_values[2] = 10000;


	expected_answer[0] = 100; expected_answer[1] = 10000; expected_answer[2] = 1000; expected_answer[3] = 1000;

	gather_values_by_deme(indices, demes, deme_specific_values, device_answer);
	thrust::host_vector<int> host_answer = device_answer;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

