#include <util/which_function.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class WhichGreaterThanIntegerSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(WhichGreaterThanIntegerSuccess, which_greater_than) {
	thrust::device_vector<int> stencil(5);
	int value = 10;
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(2);

	stencil[0] = 3; stencil[1] = 27; stencil[2] = 3; stencil[3] = 14812; stencil[4] = 3;
	which_greater_than(stencil, device_answer, value);
	thrust::host_vector<int> host_answer = device_answer;

	expected_answer[0] = 1; expected_answer[1] = 3; 

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

