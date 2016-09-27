#include <util/which_function.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class WhichGreaterThanSubsetIntegerSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(WhichGreaterThanSubsetIntegerSuccess, which_greater_than) {
	thrust::device_vector<int> stencil(5);
	int value = 10;
	int elements_evaluated = 2;
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(1);

	stencil[0] = 3; stencil[1] = 27; stencil[2] = 3; stencil[3] = 14812; stencil[4] = 3;
	which_greater_than(stencil, device_answer, value, elements_evaluated);
	thrust::host_vector<int> host_answer = device_answer;

	expected_answer[0] = 1; 

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

