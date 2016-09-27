#include <util/one_dim_two_dim.h>
#include <thrust/sequence.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class OneDimTwoDimSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(OneDimTwoDimSuccess, one_dim_two_dim) {
	thrust::device_vector<int> vector1_vals(2);
	thrust::device_vector<int> vector2_vals(3);
	thrust::device_vector<int> new_vector(6);
	thrust::device_vector<int> device_v1_answer(6);
	thrust::device_vector<int> device_v2_answer(6);

	thrust::host_vector<int> expected_v1_answer(6);
	thrust::host_vector<int> expected_v2_answer(6);

	vector1_vals[0] = 0; vector1_vals[1] = 2;
	vector2_vals[0] = 1; vector2_vals[1] = 3; vector2_vals[2] = 5;

	thrust::sequence(new_vector.begin(), new_vector.end());
	one_dim_two_dim(vector1_vals, vector2_vals, new_vector, device_v1_answer, device_v2_answer);
	thrust::host_vector<int> host_answer_1 = device_v1_answer;
	thrust::host_vector<int> host_answer_2 = device_v2_answer;

	expected_v1_answer[0] = 0; expected_v1_answer[1] = 0; expected_v1_answer[2] = 0; expected_v1_answer[3] = 2; expected_v1_answer[4] = 2; expected_v1_answer[5] = 2;

	expected_v2_answer[0] = 1; expected_v2_answer[1] = 3; expected_v2_answer[2] = 5; expected_v2_answer[3] = 1; expected_v2_answer[4] = 3; expected_v2_answer[5] = 5;

	EXPECT_THAT(host_answer_1, ::testing::ContainerEq(expected_v1_answer));
	EXPECT_THAT(host_answer_2, ::testing::ContainerEq(expected_v2_answer));
	}

