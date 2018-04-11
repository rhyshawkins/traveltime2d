
#include <cmath>
#include <gtest/gtest.h>

#include "subfield.hpp"

TEST(SubfieldSuite, HorizontalRamp)
{
  constexpr int WIDTH = 16;
  constexpr int HEIGHT = 16;
  constexpr int SIZE = WIDTH * HEIGHT;
  constexpr int SUBWIDTH = 7;
  constexpr int SUBHEIGHT = 7;

  double image[SIZE];

  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      image[j * WIDTH + i] = (double)i/(double)(WIDTH - 1);
      // printf("%5.3f ", image[j * WIDTH + i]);
    }
    // printf("\n");
  }

  Subfield<double> subfield(WIDTH, HEIGHT, 0, 0, SUBWIDTH, SUBHEIGHT);

  subfield.update(image);

  for (int j = 0; j < SUBHEIGHT; j ++) {
    for (int i = 0; i < SUBWIDTH; i ++) {
    
      EXPECT_DOUBLE_EQ(subfield.field[j * SUBWIDTH + i], (double)i/(double)(WIDTH + WIDTH - 2));
      // printf("%6.3f ", subfield.field[j * SUBWIDTH + i]);
    }
    // printf("\n");
  }

}

TEST(SubfieldSuite, VerticalRamp)
{
  constexpr int WIDTH = 16;
  constexpr int HEIGHT = 16;
  constexpr int SIZE = WIDTH * HEIGHT;
  constexpr int SUBWIDTH = 7;
  constexpr int SUBHEIGHT = 7;

  double image[SIZE];

  for (int j = 0; j < HEIGHT; j ++) {
    for (int i = 0; i < WIDTH; i ++) {

      image[j * WIDTH + i] = (double)j/(double)(HEIGHT - 1);
      // printf("%5.3f ", image[j * WIDTH + i]);
    }
    // printf("\n");
  }

  Subfield<double> subfield(WIDTH, HEIGHT, 0, 0, SUBWIDTH, SUBHEIGHT);

  subfield.update(image);

  for (int j = 0; j < SUBHEIGHT; j ++) {
    for (int i = 0; i < SUBWIDTH; i ++) {
    
      EXPECT_DOUBLE_EQ(subfield.field[j * SUBWIDTH + i], (double)j/(double)(HEIGHT + HEIGHT - 2));
      // printf("%6.3f ", subfield.field[j * SUBWIDTH + i]);
    }
    // printf("\n");
  }

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
