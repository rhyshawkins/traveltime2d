
#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

constexpr double LONMIN = -25.6;
constexpr double LONMAX = -12.4;
constexpr double LATMIN = 61.9;
constexpr double LATMAX = 67.9;

TEST(TravelTimeFieldSuite, LinearRamp65Refinement0)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 2.5 + nlon;
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.25 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					0);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					0);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Traveltime %10.6f %10.6f\n",
	 traveltimeAB,
	 traveltimeBA);
}

TEST(TravelTimeFieldSuite, Sin2Pi65Refinement0)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 3.0 + 0.5*sin(2.0 * M_PI * nlon);
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.4 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					0);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					0);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Sin2Pi Traveltime %10.6f %10.6f\n",
	 traveltimeAB,
	 traveltimeBA);
}

TEST(TravelTimeFieldSuite, Sin32Pi65Refinement0)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int REFINE = 1;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 3.0 + 0.5*sin(32.0 * M_PI * nlon);
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.4 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					REFINE);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					REFINE);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Sin32Pi Traveltime %10.6f %10.6f (%10.6f)\n",
	 traveltimeAB,
	 traveltimeBA,
	 (traveltimeAB +
	  traveltimeBA)/2.0);
}

TEST(TravelTimeFieldSuite, Sin32Pi76Refinement0)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int REFINE = 0;
  
  constexpr int SIZE = WIDTH * HEIGHT;
  
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 3.0 + 0.5*sin(32.0 * M_PI * nlon);
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.4 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					REFINE);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					REFINE);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Sin32Pi Traveltime %10.6f %10.6f\n",
	 traveltimeAB,
	 traveltimeBA);
}

TEST(TravelTimeFieldSuite, Sin32Pi87Refinement0)
{
  constexpr int WIDTH = 256;
  constexpr int HEIGHT = 128;
  constexpr int REFINE = 0;
  
  constexpr int SIZE = WIDTH * HEIGHT;
  
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 3.0 + 0.5*sin(32.0 * M_PI * nlon);
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.4 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					REFINE);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					REFINE);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Sin32Pi Traveltime %10.6f %10.6f\n",
	 traveltimeAB,
	 traveltimeBA);
}

TEST(TravelTimeFieldSuite, Sin8Pi98Refinement0)
{
  constexpr int WIDTH = 512;
  constexpr int HEIGHT = 256;
  constexpr int REFINE = 0;
  
  constexpr int SIZE = WIDTH * HEIGHT;
  
  
  double image[SIZE];
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      double nlon = ((double)i + 0.5)/(double)WIDTH;

      image[j*WIDTH + i] = 3.0 + 0.5*sin(32.0 * M_PI * nlon);
    }
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double dlon = 0.4 * (LONMAX - LONMIN);
  double clat = 0.5 * (LATMAX + LATMIN);
  LonLat<> A = LonLat<>(LONMIN + dlon, clat);
  LonLat<> B = LonLat<>(LONMAX - dlon, clat);

  TravelTimeField<LonLat<>, double> ttA(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					A,
					velocity,
					REFINE);

  ttA.construct_traveltime_field();

  double traveltimeAB = ttA.get_traveltime(B);

  
  TravelTimeField<LonLat<>, double> ttB(LonLat<>(LONMIN, LATMIN),
					LonLat<>(LONMAX, LATMAX),
					B,
					velocity,
					REFINE);

  ttB.construct_traveltime_field();

  double traveltimeBA = ttB.get_traveltime(A);

  
  printf("Sin32Pi Traveltime %10.6f %10.6f\n",
	 traveltimeAB,
	 traveltimeBA);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
      

						     

    
  
  



  
