#include <util/amplify.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class AmplifySequenceTest : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(AmplifySequenceTest, amplify_sequence) {
	thrust::device_vector<int> amplify_counts(3);
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(11);

	amplify_counts[0] = 5; amplify_counts[1] = 2; amplify_counts[2] = 4;
	int number_of_elements_in_sequence = 3;

	expected_answer[0] = 0; expected_answer[1] = 0; expected_answer[2] = 0; expected_answer[3] = 0; expected_answer[4] = 0; expected_answer[5] = 1; expected_answer[6] = 1; expected_answer[7] = 2; expected_answer[8] = 2; expected_answer[9] = 2; expected_answer[10] = 2;

	amplify_sequence(amplify_counts, number_of_elements_in_sequence, device_answer);
	thrust::host_vector<int> host_answer = device_answer;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

