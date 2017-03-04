#include <util/remove_duplicate_pairs.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class RemoveDuplicatePairsTest : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(RemoveDuplicatePairsTest, remove_duplicate_pairs) {
	thrust::device_vector<int> vector_A(5);
	thrust::device_vector<int> vector_B(5);

	vector_A[0] = 0; vector_A[1] = 0; vector_A[2] = 1; vector_A[3] = 1; vector_A[4] = 1; 
	vector_B[0] = 0; vector_B[1] = 7; vector_B[2] = 8; vector_B[3] = 9; vector_B[4] = 8; 
	
	thrust::host_vector<int> expected_answer_A(4);
	expected_answer_A[0] = 0;
	expected_answer_A[1] = 0;
	expected_answer_A[2] = 1;
	expected_answer_A[3] = 1;

	thrust::host_vector<int> expected_answer_B(4);
	expected_answer_B[0] = 0;
	expected_answer_B[1] = 7;
	expected_answer_B[2] = 8;
	expected_answer_B[3] = 9;

	remove_duplicate_pairs(vector_A, vector_B);
	thrust::host_vector<int> host_answer_A = vector_A;
	thrust::host_vector<int> host_answer_B = vector_B;

	EXPECT_THAT(host_answer_A, ::testing::ContainerEq(expected_answer_A));
	EXPECT_THAT(host_answer_B, ::testing::ContainerEq(expected_answer_B));
	}

