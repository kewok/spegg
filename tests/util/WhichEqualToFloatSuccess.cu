#include <util/which_function.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// The fixture for testing class Foo.
class WhichEqualToFloatSuccess : public ::testing::Test {
	protected:
		virtual void setup() {
		}
		
		virtual void TearDown() {
		}
};

TEST_F(WhichEqualToFloatSuccess, which_equal_to) {
	thrust::device_vector<float> stencil(5);
	float value = 3.14;
	thrust::device_vector<int> device_answer;
	thrust::host_vector<int> expected_answer(3);

	stencil[0] = 3.14; stencil[1] = 2.5; stencil[2] = 3.14; stencil[3] = 14812.0487; stencil[4] = 3.14;
	which_equal_to(stencil, device_answer, value);
	thrust::host_vector<int> host_answer = device_answer;

	expected_answer[0] = 0; expected_answer[1] = 2; expected_answer[2] = 4;

	EXPECT_THAT(host_answer, ::testing::ContainerEq(expected_answer));
	}

