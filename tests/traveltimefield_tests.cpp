

#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

TEST(TravelTimeFieldSuite, SmallHomogenous) {

  double image[4];
  VelocityField<double> velocity(image, 2, 2);

  //
  // Initialize constant image
  //
  image[0] = 1.0;
  image[1] = 1.0;
  image[2] = 1.0;
  image[3] = 1.0;

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  0);

  tt.construct_traveltime_field();

  for (int i = 0; i < 4; i ++) {
    for (int j = 0; j < 4; j ++) {

      int o = i * 4 + j;
      EXPECT_EQ(tt.nodes[o].state, TravelTimeNode<double>::STATE_FIXED);

      EXPECT_EQ(tt.nodes[o].nx, (1.5 * j)/3.0 - 0.25);
      EXPECT_EQ(tt.nodes[o].ny, (1.5 * i)/3.0 - 0.25);
    }
  }

  constexpr double A = 1.0/sqrt(2.0);
  constexpr double B = A + 1.0;
  constexpr double C = B + A;
  constexpr double D = ((A + B)/2.0 + (B + C)/2.0)/2.0;

  //
  // Top Row
  //
  EXPECT_DOUBLE_EQ(tt.nodes[0].T, C);
  EXPECT_DOUBLE_EQ(tt.nodes[1].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[2].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[3].T, C);

  //
  // Middle Rows
  //
  EXPECT_DOUBLE_EQ(tt.nodes[4].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[5].T, A);
  EXPECT_DOUBLE_EQ(tt.nodes[6].T, A);
  EXPECT_DOUBLE_EQ(tt.nodes[7].T, B);
  
  EXPECT_DOUBLE_EQ(tt.nodes[8].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[9].T, A);
  EXPECT_DOUBLE_EQ(tt.nodes[10].T, A);
  EXPECT_DOUBLE_EQ(tt.nodes[11].T, B);

  //
  // Bottom Row
  //
  EXPECT_DOUBLE_EQ(tt.nodes[12].T, C);
  EXPECT_DOUBLE_EQ(tt.nodes[13].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[14].T, B);
  EXPECT_DOUBLE_EQ(tt.nodes[15].T, C);

  //
  // Interpolated travel times
  //
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-0.1, -0.1)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-0.1,  0.1)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 0.1,  0.1)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 0.1, -0.1)), A);

  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-0.5, -0.5)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-0.5,  0.5)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 0.5,  0.5)), A);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 0.5, -0.5)), A);

  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-1.0, -1.0)), D);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm(-1.0,  1.0)), D);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 1.0,  1.0)), D);
  EXPECT_DOUBLE_EQ(tt.get_traveltime(CartesianKm( 1.0, -1.0)), D);
}

TEST(TravelTimeFieldSuite, LonLatHomogenous)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  
  double *image = new double[WIDTH * HEIGHT];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      image[j * WIDTH + i] = 3.0;
    }
  }
  
  TravelTimeField<LonLat<>, double> tt(LonLat<>(-10.0, -10.0),
				       LonLat<>(10.0, 10.0),
				       LonLat<>(0.0, 0.0),
				       velocity,
				       0);

  // for (size_t i = 0; i < tt.width; i ++) {
  //   printf("%d %f\n", (int)i, tt.nodes[i].nx);
  // }

  tt.construct_traveltime_field();

  if (!tt.save_traveltime_field("LonLatHomogenous.txt")) {
    fprintf(stderr, "error: failed to save traveltime field\n");
  }
  
  printf("%f %f\n",
	 LonLat<>::distance_km(LonLat<>(0.0, -10.0), LonLat<>(0.0, 0.0)),
	 LonLat<>::distance_km(LonLat<>(0.0, 10.0), LonLat<>(0.0, 0.0)));
			       
  printf("%f\n", tt.get_traveltime(LonLat<>(0.0, -10.0)));
  printf("%f\n", tt.get_traveltime(LonLat<>(0.0, 10.0)));
  
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, -7.5)));
  

}

TEST(TravelTimeFieldSuite, LonLatHomogenousHighLat)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  
  double *image = new double[WIDTH * HEIGHT];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      image[j * WIDTH + i] = 3.0;
    }
  }
  
  TravelTimeField<LonLat<>, double> tt(LonLat<>(-10.0, 50.0),
				       LonLat<>(10.0, 70.0),
				       LonLat<>(0.0, 60.0),
				       velocity,
				       0);

  // for (size_t i = 0; i < tt.width; i ++) {
  //   printf("%d %f\n", (int)i, tt.nodes[i].nx);
  // }

  tt.construct_traveltime_field();

  if (!tt.save_traveltime_field("LonLatHomogenousHighLat.txt")) {
    fprintf(stderr, "error: failed to save traveltime field\n");
  }
  
  printf(" *** High Lat ***\n");
  printf("Distance bottom to centre %f top to centre %f\n",
	 LonLat<>::distance_km(LonLat<>(0.0, 50.0), LonLat<>(0.0, 60.0)),
	 LonLat<>::distance_km(LonLat<>(0.0, 70.0), LonLat<>(0.0, 60.0)));
			       
  printf("Time to bottom %f\n", tt.get_traveltime(LonLat<>(0.0, 50.0)));
  printf("Time to top %f\n", tt.get_traveltime(LonLat<>(0.0, 70.0)));
  
  printf("Distance Left to centre %f Right to centre %f\n",
	 LonLat<>::distance_km(LonLat<>(10.0, 60.0), LonLat<>(0.0, 60.0)),
	 LonLat<>::distance_km(LonLat<>(-10.0, 60.0), LonLat<>(0.0, 60.0)));
			       
  printf("Time to left %f\n", tt.get_traveltime(LonLat<>(10.0, 60.0)));
  printf("Time to right %f\n", tt.get_traveltime(LonLat<>(-10.0, 60.0)));

  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, -7.5)));
  

}

TEST(TravelTimeFieldSuite, LonLatHomogenousLowLat)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  
  double *image = new double[WIDTH * HEIGHT];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      image[j * WIDTH + i] = 3.0;
    }
  }
  
  TravelTimeField<LonLat<>, double> tt(LonLat<>(-10.0, -70.0),
				       LonLat<>(10.0, -50.0),
				       LonLat<>(0.0, -60.0),
				       velocity,
				       0);

  // for (size_t i = 0; i < tt.width; i ++) {
  //   printf("%d %f\n", (int)i, tt.nodes[i].nx);
  // }

  tt.construct_traveltime_field();

  if (!tt.save_traveltime_field("LonLatHomogenousLowLat.txt")) {
    fprintf(stderr, "error: failed to save traveltime field\n");
  }

  printf(" *** Low Lat ***\n");
  printf("Distance bottom to centre %f top to centre %f\n",
	 LonLat<>::distance_km(LonLat<>(0.0, -70.0), LonLat<>(0.0, -60.0)),
	 LonLat<>::distance_km(LonLat<>(0.0, -50.0), LonLat<>(0.0, -60.0)));
			       
  printf("Time to bottom %f\n", tt.get_traveltime(LonLat<>(0.0, -70.0)));
  printf("Time to top %f\n", tt.get_traveltime(LonLat<>(0.0, -50.0)));
  
  printf("Distance Left to centre %f Right to centre %f\n",
	 LonLat<>::distance_km(LonLat<>(10.0, -60.0), LonLat<>(0.0, -60.0)),
	 LonLat<>::distance_km(LonLat<>(-10.0, -60.0), LonLat<>(0.0, -60.0)));
			       
  printf("Time to left %f\n", tt.get_traveltime(LonLat<>(10.0, -60.0)));
  printf("Time to right %f\n", tt.get_traveltime(LonLat<>(-10.0, -60.0)));

  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(-7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, 0.0)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(7.5, -7.5)));
  // printf("%f\n", tt.get_traveltime(LonLat<>(0.0, -7.5)));
  

}

TEST(TravelTimeFieldSuite, MediumHomogenousRefinement1)
{
  constexpr int WIDTH = 32;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 1.0;
  }

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  1);

  

  tt.construct_traveltime_field();

  for (int i = 0; i < WIDTH; i ++) {
    for (int j = 0; j < HEIGHT; j ++) {

      int o = i * WIDTH + j;
      EXPECT_EQ(tt.nodes[o].state, TravelTimeNode<double>::STATE_FIXED);
    }
  }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
