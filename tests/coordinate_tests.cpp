

#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"

TEST(CoordinateSuite, Spherical) {

  LonLat<> a(0.0, 0.0);
  LonLat<> b(90.0, 0.0);
  
  double dist = LonLat<>::distance_km(a, b);

  EXPECT_DOUBLE_EQ(dist, 0.5 * 6371.0 * M_PI);
  
}

TEST(CoordinateSuite, Cartesian) {

  CartesianKm a(0.0, 0.0);
  CartesianKm b(1.0, 1.0);
  
  double dist = CartesianKm::distance_km(a, b);

  EXPECT_DOUBLE_EQ(dist, sqrt(2.0));

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
