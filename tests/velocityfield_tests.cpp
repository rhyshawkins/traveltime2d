

#include <cmath>
#include <gtest/gtest.h>

#include "velocityfield.hpp"

TEST(VelocityFieldSuite, Lerp) {

  double image[4];
  VelocityField<double> field(image, 2, 2);

  //
  // Initialize constant image
  //
  image[0] = 1.0;
  image[1] = 1.0;
  image[2] = 1.0;
  image[3] = 1.0;

  constexpr int N = 10;

  for (int i = 0; i < N; i ++) {

    double nx = (double)i/(double)(N - 1);
    
    for (int j = 0; j < N; j ++) {

      double ny = (double)j/(double)(N - 1);

      EXPECT_DOUBLE_EQ(field.lerp(nx, ny), 1.0);

    }
  }

  //
  // Initialize horizontal gradient
  //
  image[0] = 1.0;
  image[1] = 2.0;
  image[2] = 1.0;
  image[3] = 2.0;

  for (int i = 0; i < N; i ++) {

    double nx = (double)i/(double)(N - 1);

    double expected;
    if (nx < 0.25) {
      expected = 1.0;
    } else if (nx > 0.75) {
      expected = 2.0;
    } else {
      expected = 1.0 + (nx - 0.25)/0.5;
    }
    
    for (int j = 0; j < N; j ++) {

      double ny = (double)j/(double)(N - 1);

      // printf("%10.6f %10.6f %10.6f %10.6f\n", nx, ny, field.lerp(nx, ny), expected);

      EXPECT_DOUBLE_EQ(field.lerp(nx, ny), expected);

    }
  }
      
}

TEST(VelocityFieldSuite, LerpLarge) {

  double image[128*64];
  VelocityField<double> field(image, 128, 64);

  //
  // Initialize constant image
  //
  for (int j = 0; j < 64; j ++) {
    for (int i = 0; i < 128; i ++) {
      image[j * 128 + i] = 3.0;
    }
  }

  constexpr int M = 20;
  constexpr int N = 10;

  for (int i = 0; i < M; i ++) {

    double nx = (double)i/(double)(M - 1);
    
    for (int j = 0; j < N; j ++) {

      double ny = (double)j/(double)(N - 1);

      EXPECT_DOUBLE_EQ(field.lerp(nx, ny), 3.0);

    }
  }
  
  //
  // Initialize perturbation
  //
  for (int j = 32; j < 64; j ++) {
    for (int i = 32; i < 64; i ++) {
      image[j * 128 + i] = 2.5;
    }
  }

  for (int j = 0; j < N; j ++) {
    
    double ny = (double)j/(double)(N - 1);
    for (int i = 0; i < M; i ++) {
      
      double nx = (double)i/(double)(M - 1);

      if (j >= N/2 && i >= M/4 && i < M/2) {
	EXPECT_DOUBLE_EQ(field.lerp(nx, ny), 2.5);
      } else {
	EXPECT_DOUBLE_EQ(field.lerp(nx, ny), 3.0);
      }	

    }
  }
  
  EXPECT_DOUBLE_EQ(field.lerp(0.25390625, 0.4921875), 3.0);
  EXPECT_DOUBLE_EQ(field.lerp(0.25390625, 0.5), 2.75);
  EXPECT_DOUBLE_EQ(field.lerp(0.25390625, 0.5078125), 2.5);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
