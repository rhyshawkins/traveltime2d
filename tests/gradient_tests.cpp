

#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

static double evaluate2x2uniform(int vi, double dv, int ni)
{
  double image[4];
  VelocityField<double> velocity(image, 2, 2);

  //
  // Initialize constant image
  //
  image[0] = 1.0;
  image[1] = 1.0;
  image[2] = 1.0;
  image[3] = 1.0;

  image[vi] += dv;

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  0);

  tt.construct_traveltime_field();

  return tt.nodes[ni].T;
}

static double numericaldifferentiate2x2uniform(int vi, int ni, double dv)
{
  return (evaluate2x2uniform(vi, dv/2.0, ni) -
	  evaluate2x2uniform(vi, -dv/2.0, ni))/dv;
}

static double evaluatelargeuniform(int vi, double dv, int ni, int refinement)
{
  static constexpr int WIDTH = 128;
  static constexpr int HEIGHT = 128;
  static constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 1.0;
  }

  image[vi] += dv;

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  refinement);

  tt.construct_traveltime_field();

  printf("TT: %10.6f\n", tt.nodes[ni].T);
  return tt.nodes[ni].T;
}

static double numericaldifferentiatelargeuniform(int vi, int ni, double dv, int refinement)
{
  return (evaluatelargeuniform(vi, dv/2.0, ni, refinement) -
	  evaluatelargeuniform(vi, -dv/2.0, ni, refinement))/dv;
}

TEST(TravelTimeFieldSuite, DISABLED_SmallHomogenous) {

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

  tt.nodes[0].back_project(&tt);
  tt.nodes[3].back_project(&tt);
  tt.nodes[12].back_project(&tt);
  tt.nodes[15].back_project(&tt);

  printf("Node 5 : %10.6f\n", tt.nodes[5].T);
  for (int i = 0; i < tt.nodes[5].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[5].vweights.indices[i], tt.nodes[5].vweights.weights[i]);
  }

  printf("Node 6 : %10.6f\n", tt.nodes[5].T);
  for (int i = 0; i < tt.nodes[6].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[6].vweights.indices[i], tt.nodes[6].vweights.weights[i]);
  }

  printf("Node 9 : %10.6f\n", tt.nodes[5].T);
  for (int i = 0; i < tt.nodes[9].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[9].vweights.indices[i], tt.nodes[9].vweights.weights[i]);
  }
  
  printf("Node 10: %10.6f\n", tt.nodes[5].T);
  for (int i = 0; i < tt.nodes[10].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[10].vweights.indices[i], tt.nodes[10].vweights.weights[i]);
  }

  printf("Node 4 : %10.6f\n", tt.nodes[4].T);
  for (int i = 0; i < tt.nodes[4].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[4].vweights.indices[i], tt.nodes[4].vweights.weights[i]);
  }
  printf("Node 7 : %10.6f\n", tt.nodes[7].T);
  for (int i = 0; i < tt.nodes[7].vweights.n; i ++) {
    printf("%2d %2d %10.6f\n", i, tt.nodes[7].vweights.indices[i], tt.nodes[7].vweights.weights[i]);
  }

  printf("0 7 %10.6f\n", numericaldifferentiate2x2uniform(0, 7, 1.0e-6));
  printf("1 7 %10.6f\n", numericaldifferentiate2x2uniform(1, 7, 1.0e-6));
}

TEST(TravelTimeFieldSuite, DISABLED_LargeHomogenousNumericalDiff)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 1.0;
  }
  
  // for (int j = 0; j < HEIGHT; j ++) {
  //   double v = cos(2.0 * M_PI * (double)j/(double)HEIGHT * 8.0);
  //   for (int i = 0; i < WIDTH; i ++) {
  //     double h = cos(2.0 * M_PI * (double)i/(double)WIDTH * 8.0);
  //     image[j * WIDTH + i] = 1.0 + 0.25*v*h;
  //   }
  // }

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  1);

  tt.construct_traveltime_field();

  
  tt.nodes[0].back_project(&tt);
  tt.nodes[WIDTH - 1].back_project(&tt);
  tt.nodes[(HEIGHT - 1) * WIDTH].back_project(&tt);
  tt.nodes[SIZE - 1].back_project(&tt);

  int ni = 8384;
  int vi = 8127;
  
  printf("Node %d : %15.9f %d %d\n", ni, tt.nodes[ni].T, tt.nodes[ni].vweights.n, (int)tt.nodes[ni].vweights.validate());
  double maxw = -1.0e9;
  double minw = 1.0e9;
  int minwi = 0;
  int maxwi = 0;
  for (int i = 0; i < tt.nodes[ni].vweights.n; i ++) {
    if (tt.nodes[ni].vweights.indices[i] == vi) {
      printf("%2d %2d %15.9f\n",
	     i,
	     tt.nodes[ni].vweights.indices[i],
	     tt.nodes[ni].vweights.weights[i]);
    }

    if (tt.nodes[ni].vweights.weights[i] > maxw) {
      maxw = tt.nodes[ni].vweights.weights[i];
      maxwi = tt.nodes[ni].vweights.indices[i];
    } else if (tt.nodes[ni].vweights.weights[i] < minw) {
      minw = tt.nodes[ni].vweights.weights[i];
      minwi = tt.nodes[ni].vweights.indices[i];
    }
  }

  printf("Min %2d %15.9f\n", minwi, minw);
  printf("Max %2d %15.9f\n", maxwi, maxw);
  
  printf("%5d %5d %15.9f\n", vi, ni, numericaldifferentiatelargeuniform(vi, ni, 1.0e-6, 1));
  
}

TEST(TravelTimeFieldSuite, DISABLED_LargeHomogenous)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    double v = cos(2.0 * M_PI * (double)j/(double)HEIGHT * 8.0);
    for (int i = 0; i < WIDTH; i ++) {
      double h = cos(2.0 * M_PI * (double)i/(double)WIDTH * 8.0);
      image[j * WIDTH + i] = 1.0 + 0.25*v*h;
    }
  }

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  0);

  tt.construct_traveltime_field();

  double sensitivity[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    sensitivity[i] = 0.0;
  }

  VelocityWeights<double> weights;
  double ttime = tt.get_traveltime_weighted(CartesianKm(-0.1, 0.1), weights);

  printf("Travel time: %10.6f\n", ttime);
  printf("Nweights:    %2d\n", weights.n);

  for (int i = 0; i < weights.n; i ++) {
    sensitivity[weights.indices[i]] = weights.weights[i];
  }

  FILE *fp = fopen("gradient_sensitivity.txt", "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create output file\n");
    return;
  }

  double minv = 0.0;
  int mini = 0;
  int minj = 0;
  double maxv = 0.0;
  int maxi = 0;
  int maxj = 0;
  
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      fprintf(fp, "%16.9e ", sensitivity[j * WIDTH + i]);

      if (sensitivity[j * WIDTH + i] < minv) {
	minv = sensitivity[j * WIDTH + i];
	mini = i;
	minj = j;
      } else if (sensitivity[j * WIDTH + i] > maxv) {
	maxv = sensitivity[j * WIDTH + i];
	maxi = i;
	maxj = j;
      }
	
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  printf("Min: %10.6f (%2d, %2d) %d\n", minv, mini, minj, minj * WIDTH + mini);
  printf("Max: %10.6f (%2d, %2d) %d\n", maxv, maxi, maxj, maxj * WIDTH + maxi);

  fp = fopen("gradient_image.txt", "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create output file\n");
    return;
  }
  
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      fprintf(fp, "%16.9e ", image[j * WIDTH + i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}


TEST(TravelTimeFieldSuite, LargeHomogenousRefined)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    double v = cos(2.0 * M_PI * (double)j/(double)HEIGHT * 8.0);
    for (int i = 0; i < WIDTH; i ++) {
      double h = cos(2.0 * M_PI * (double)i/(double)WIDTH * 8.0);
      image[j * WIDTH + i] = 1.0 + 0.25*v*h;
    }
  }

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(0.0, 0.0),
					  velocity,
					  1);

  tt.construct_traveltime_field();

  double sensitivity[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    sensitivity[i] = 0.0;
  }

  VelocityWeights<double> weights;
  double ttime = tt.get_traveltime_weighted(CartesianKm(-0.5, 0.5), weights);

  printf("Travel time: %10.6f\n", ttime);
  printf("Nweights:    %2d\n", weights.n);

  for (int i = 0; i < weights.n; i ++) {
    sensitivity[weights.indices[i]] = weights.weights[i];
  }

  FILE *fp = fopen("gradient_sensitivity_refined.txt", "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create output file\n");
    return;
  }
  
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      fprintf(fp, "%16.9e ", sensitivity[j * WIDTH + i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

}

TEST(TravelTimeFieldSuite, DISABLED_LargeHomogenousRefinedReverse)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize constant image
  //
  for (int j = 0; j < HEIGHT; j ++) {
    double v = cos(2.0 * M_PI * (double)j/(double)HEIGHT * 8.0);
    for (int i = 0; i < WIDTH; i ++) {
      double h = cos(2.0 * M_PI * (double)i/(double)WIDTH * 8.0);
      image[j * WIDTH + i] = 1.0 + 0.25*v*h;
    }
  }

  TravelTimeField<CartesianKm, double> tt(CartesianKm(-1.0, -1.0),
					  CartesianKm(1.0, 1.0),
					  CartesianKm(-0.5, 0.5),
					  velocity,
					  1);

  tt.construct_traveltime_field();

  double sensitivity[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    sensitivity[i] = 0.0;
  }

  VelocityWeights<double> weights;
  double ttime = tt.get_traveltime_weighted(CartesianKm(0.0, 0.0), weights);

  printf("Travel time: %10.6f\n", ttime);
  printf("Nweights:    %2d\n", weights.n);

  for (int i = 0; i < weights.n; i ++) {
    sensitivity[weights.indices[i]] = weights.weights[i];
  }

  FILE *fp = fopen("gradient_sensitivity_refined_reverse.txt", "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create output file\n");
    return;
  }
  
  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {
      fprintf(fp, "%16.9e ", sensitivity[j * WIDTH + i]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

