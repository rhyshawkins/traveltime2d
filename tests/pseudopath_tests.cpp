

#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

int save_image(const char *filename, double *img, int width, int height)
{
  FILE *fp;
  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%15.9f ", img[j * width + i]);

    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

TEST(TravelTimeFieldSuite, Gaussian) {

  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 128;
  
  double image[WIDTH * HEIGHT];
  double TTa[(WIDTH + 2) * (HEIGHT + 2)];
  double TTb[(WIDTH + 2) * (HEIGHT + 2)];
    
  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  //
  // Initialize slow Gaussian hill in centre
  //
  for (int j = 0; j < HEIGHT; j ++) {

    double lat = -1.0 + 2.0 * (double)j/(double)(HEIGHT - 1);
    
    for (int i = 0; i < WIDTH; i ++) {

      double lon = -1.0 + 2.0 * (double)i/(double)(WIDTH - 1);

      image[j * WIDTH + i] = 3.0 - 2.5*exp(-(lon*lon + lat*lat)/(2.0 * 0.25*0.25));
    }
  }

  CartesianKm A(-0.5, 0.0);
  CartesianKm B(0.5, 0.0);
  
  TravelTimeField<CartesianKm, double> ttAB(CartesianKm(-1.0, -1.0),
					    CartesianKm(1.0, 1.0),
					    A,
					    velocity,
					    0);

  ttAB.construct_traveltime_field();

  TravelTimeField<CartesianKm, double> ttBA(CartesianKm(-1.0, -1.0),
					    CartesianKm(1.0, 1.0),
					    B,
					    velocity,
					    0);

  ttBA.construct_traveltime_field();

  
  printf("%f\n", ttAB.get_traveltime(B));
  printf("%f\n", ttBA.get_traveltime(A));

  ttAB.get_traveltime_image(TTa);
  save_image("pseudopath_image_A.txt", TTa, WIDTH + 2, HEIGHT + 2);
  ttBA.get_traveltime_image(TTb);
  save_image("pseudopath_image_B.txt", TTb, WIDTH + 2, HEIGHT + 2);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
