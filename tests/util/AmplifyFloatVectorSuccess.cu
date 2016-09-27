#include <util/amplify.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class AmplifyFloatVectorTest : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(AmplifyFloatVectorTest, amplify_float) {
	thrust::device_vector<float> values_to_amplify(3);
	thrust::device_vector<int> amplify_counts(3);
	thrust::device_vector<float> device_answer;
	thrust::host_vector<float> expected_answer(11);

	values_to_amplify[0] = 3.14; values_to_amplify[1] = 2.73; values_to_amplify[2] = 37.58;
	
	amplify_counts[0] = 5; amplify_counts[1] = 2; amplify_counts[2] = 4;

	expected_answer[0] = 3.14; expected_answer[1] = 3.14; expected_answer[2] = 3.14; expected_answer[3] = 3.14; expected_answer[4] = 3.14; expected_answer[5] = 2.73; expected_answer[6] = 2.73; expected_answer[7] = 37.58; expected_answer[8] = 37.58; expected_answer[9] = 37.58; expected_answer[10] = 37.58;

	amplify_float(values_to_amplify, amplify_counts, device_answer);
	thrust::host_vector<float> host_answer = device_answer;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

