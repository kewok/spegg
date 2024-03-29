#include <util/reduce_by_key_with_zeroes.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class ReduceByKeywZerosIntegerTest : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(ReduceByKeywZerosIntegerTest, reduce_by_key_with_zeros) {
	thrust::device_vector<int> keys(5);
	thrust::device_vector<int> values(5);
	thrust::device_vector<int> device_answer;
	int number_of_values = 5;
	int number_of_possible_keys = 5;

	keys[0]=1; keys[1] = 3; keys[2] = 4; keys[3] = 4; keys[4] = 4;
	values[0] = 5; values[1] = 6; values[2] =  7; values[3] = 8; values[4] = 9;
	
	thrust::host_vector<int> expected_answer(5);
	expected_answer[0] = 0;
	expected_answer[1] = 5;
	expected_answer[2] = 0;
	expected_answer[3] = 6;
	expected_answer[4] = 24;

	reduce_by_key_with_zeros(keys, values, device_answer, number_of_values, number_of_possible_keys);
	thrust::host_vector<int> host_answer = device_answer;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

