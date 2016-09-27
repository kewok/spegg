#include <util/amplify.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class AmplifyIntegerVectorTest : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(AmplifyIntegerVectorTest, amplify) {
	thrust::device_vector<int> values_to_amplify(3);
	thrust::device_vector<int> amplify_counts(3);
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(11);

	values_to_amplify[0] = 1; values_to_amplify[1] = 2; values_to_amplify[2] = 3;
	
	amplify_counts[0] = 5; amplify_counts[1] = 2; amplify_counts[2] = 4;

	expected_answer[0] = 1; expected_answer[1] = 1; expected_answer[2] = 1; expected_answer[3] = 1; expected_answer[4] = 1; expected_answer[5] = 2; expected_answer[6] = 2; expected_answer[7] = 3; expected_answer[8] = 3; expected_answer[9] = 3; expected_answer[10] = 3;

	amplify(values_to_amplify, amplify_counts, device_answer);
	thrust::host_vector<int> host_answer = device_answer;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

