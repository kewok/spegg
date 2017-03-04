#include <util/which_function.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class WhichLessThanFloatSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(WhichLessThanFloatSuccess, which_less_than) {
	thrust::device_vector<float> stencil(5);
	float value = 10.547;
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(3);

	stencil[0] = 10.54; stencil[1] = 27.24; stencil[2] = 32.24; stencil[3] = 1.4812; stencil[4] = 0.3;
	which_less_than(stencil, device_answer, value);
	thrust::host_vector<int> host_answer = device_answer;

	expected_answer[0] = 0; expected_answer[1] = 3; expected_answer[2] = 4;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

