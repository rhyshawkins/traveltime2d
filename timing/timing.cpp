
#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

extern "C" {
#include "tracking.h"
};

static char short_options[] = "x:y:z:h";

static struct option long_options[] = {
  {"width", required_argument, 0, 'x'},
  {"height", required_argument, 0, 'y'},
  {"depth", required_argument, 0, 'z'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  int width = 128;
  int height = 64;
  int size;
  int iterations = 100;

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'x':
      width = atoi(optarg);
      break;

    case 'y':
      height = atoi(optarg);
      break;

    case 'z':
      iterations = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }
  
  double *image;

  size = width * height;
  image = new double[size];

  for (int i = 0; i < size; i ++) {
    image[i] = 1.0;
  }
  
  VelocityField<double> *velocity = new VelocityField<double>(image, width, height);
  
  TravelTimeField<LonLat<>, double> tt(LonLat<>(-10.0, -10.0),
				       LonLat<>(10.0, 10.0),
				       LonLat<>(0.0, 0.0),
				       *velocity,
				       0);
  tracking_t makefield;

  tracking_init(&makefield);
    
  for (int i = 0; i < iterations; i ++) {
    tracking_start(&makefield);
    tt.construct_traveltime_field();
    tracking_end(&makefield);
  }

  printf("%d x %d x %d\n", width, height, iterations);

  double mean = tracking_mean(&makefield);
  printf("us/iteration  : %10.6f us\n", mean);

  mean = mean/1.0e6;
  printf("1k forward    :  %10.6f hr\n", (mean * iterations) * 1000.0/3600.0);
  printf("1M forward    : %10.6f hr\n", (mean * iterations) * 1000000.0/3600.0);
  printf("1k forward/16 : %10.6f hr\n", ((mean * iterations) * 1000.0/3600.0)/16.0);
  printf("1M forward/16 : %10.6f hr\n", ((mean * iterations) * 1000000.0/3600.0)/16.0);
  printf("1k forward/12 : %10.6f hr\n", ((mean * iterations) * 1000.0/3600.0)/12.0);
  printf("1M forward/12 : %10.6f hr\n", ((mean * iterations) * 1000000.0/3600.0)/12.0);
  
  return 0;
}


static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of\n"
	  "\n"
	  " -x|--width <int>         Pixels in x direction\n"
	  " -y|--height <int>        Pixels in y direction\n"
	  " -z|--depth <int>         No. iterations\n"
	  "\n"
	  " -h|--help                Show usage\n"
	  "", pname);
}
