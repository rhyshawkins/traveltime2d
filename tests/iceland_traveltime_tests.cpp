
#include <cmath>
#include <gtest/gtest.h>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimefield.hpp"

constexpr double LONMIN = -25.6;
constexpr double LONMAX = -12.4;
constexpr double LATMIN = 61.9;
constexpr double LATMAX = 67.9;

std::vector<LonLat<>> STATIONS = {
  LonLat<>(-21.903950000, 64.127670000),
  LonLat<>(-21.326800000, 64.747400000),
  LonLat<>(-18.099773000, 65.685928000),
  LonLat<>(-19.599640000, 65.670700000),
  LonLat<>(-20.721770000, 65.422218000),
  LonLat<>(-22.428350000, 65.926750000),
  LonLat<>(-14.501803000, 65.256226000),
  LonLat<>(-13.753478000, 65.540825000),
  LonLat<>(-15.169086000, 66.121422000),
  LonLat<>(-18.256500000, 65.302811000),
  LonLat<>(-24.161409000, 65.609856000),
  LonLat<>(-14.090622000, 64.810371000),
  LonLat<>(-15.308651000, 65.165932000),
  LonLat<>(-23.486561000, 65.873932000),
  LonLat<>(-22.232260000, 64.745897000),
  LonLat<>(-22.510008000, 65.598396000),
  LonLat<>(-18.331696000, 65.027626000),
  LonLat<>(-19.673590000, 65.225365000),
  LonLat<>(-15.353522000, 64.886208000),
  LonLat<>(-16.650196000, 65.054375000),
  LonLat<>(-18.130611000, 63.769817000),
  LonLat<>(-17.266495000, 64.406776000),
  LonLat<>(-15.139223000, 64.287781000),
  LonLat<>(-16.641010000, 63.876549000),
  LonLat<>(-22.423315000, 65.180504000),
  LonLat<>(-21.096394000, 65.109688000),
  LonLat<>(-21.677841000, 65.705025000),
  LonLat<>(-14.873623000, 65.728622000),
  LonLat<>(-21.167883000, 64.494139000),
  LonLat<>(-19.484184000, 64.532478000),
  LonLat<>(-23.852322000, 64.907349000)
};

std::vector<CartesianKm> CARTESIAN_STATIONS = {
  CartesianKm(-21.903950000, 64.127670000),
  CartesianKm(-21.326800000, 64.747400000),
  CartesianKm(-18.099773000, 65.685928000),
  CartesianKm(-19.599640000, 65.670700000),
  CartesianKm(-20.721770000, 65.422218000),
  CartesianKm(-22.428350000, 65.926750000),
  CartesianKm(-14.501803000, 65.256226000),
  CartesianKm(-13.753478000, 65.540825000),
  CartesianKm(-15.169086000, 66.121422000),
  CartesianKm(-18.256500000, 65.302811000),
  CartesianKm(-24.161409000, 65.609856000),
  CartesianKm(-14.090622000, 64.810371000),
  CartesianKm(-15.308651000, 65.165932000),
  CartesianKm(-23.486561000, 65.873932000),
  CartesianKm(-22.232260000, 64.745897000),
  CartesianKm(-22.510008000, 65.598396000),
  CartesianKm(-18.331696000, 65.027626000),
  CartesianKm(-19.673590000, 65.225365000),
  CartesianKm(-15.353522000, 64.886208000),
  CartesianKm(-16.650196000, 65.054375000),
  CartesianKm(-18.130611000, 63.769817000),
  CartesianKm(-17.266495000, 64.406776000),
  CartesianKm(-15.139223000, 64.287781000),
  CartesianKm(-16.641010000, 63.876549000),
  CartesianKm(-22.423315000, 65.180504000),
  CartesianKm(-21.096394000, 65.109688000),
  CartesianKm(-21.677841000, 65.705025000),
  CartesianKm(-14.873623000, 65.728622000),
  CartesianKm(-21.167883000, 64.494139000),
  CartesianKm(-19.484184000, 64.532478000),
  CartesianKm(-23.852322000, 64.907349000)
};

std::vector<std::string> STATIONNAMES = {
  "HOT30",
  "BORG ",
  "HOT13",
  "HOT12",
  "HOT11",
  "HOT10",
  "HOT17",
  "HOT16",
  "HOT15",
  "HOT14",
  "HOT08",
  "HOT19",
  "HOT18",
  "HOT09",
  "HOT02",
  "HOT07",
  "HOT26",
  "HOT27",
  "HOT24",
  "HOT25",
  "HOT22",
  "HOT23",
  "HOT20",
  "HOT21",
  "HOT04",
  "HOT05",
  "HOT06",
  "HOT29",
  "HOT01",
  "HOT28",
  "HOT03"
};

TEST(TravelTimeFieldSuite, Homogenous65Refinement0)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 0);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%2d %2d : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     (int)i, (int)j,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous65Refinement0Cubic)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 0);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime_cubic(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%2d %2d : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     (int)i, (int)j,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous65Refinement1)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 1);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%2d %2d : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     (int)i, (int)j,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous65Refinement1Cubic)
{
  constexpr int WIDTH = 64;
  constexpr int HEIGHT = 32;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 1);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime_cubic(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%2d %2d : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     (int)i, (int)j,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous76Refinement0)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 0);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%5s %5s : %10.6fkm : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     STATIONNAMES[i].c_str(), STATIONNAMES[j].c_str(),
      // 	     distkm,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous76Refinement0Cubic)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 0);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime_cubic(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%5s %5s : %10.6fkm : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     STATIONNAMES[i].c_str(), STATIONNAMES[j].c_str(),
      // 	     distkm,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous76Refinement1)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 1);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%5s %5s : %10.6fkm : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     STATIONNAMES[i].c_str(), STATIONNAMES[j].c_str(),
      // 	     distkm,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, Homogenous76Refinement1Cubic)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int SIZE = WIDTH * HEIGHT;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = 3.0;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {

    TravelTimeField<LonLat<>, double> tt(LonLat<>(LONMIN, LATMIN),
					 LonLat<>(LONMAX, LATMAX),
					 STATIONS[i],
					 velocity,
					 1);

    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime_cubic(STATIONS[j]);
      double distkm = LonLat<>::distance_km(STATIONS[i], STATIONS[j]);
      double expected_traveltime = distkm/3.0;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - 3.0);
      double relvelerr = velerr/3.0;

      // printf("%5s %5s : %10.6fkm : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      // 	     STATIONNAMES[i].c_str(), STATIONNAMES[j].c_str(),
      // 	     distkm,
      // 	     traveltime, distkm/3.0,
      // 	     distkm/traveltime, 3.0);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

TEST(TravelTimeFieldSuite, HomogenousCartesian76Refinement1Cubic)
{
  constexpr int WIDTH = 128;
  constexpr int HEIGHT = 64;
  constexpr int SIZE = WIDTH * HEIGHT;
  constexpr double TRUEVELOCITY = 1.0;
  
  double image[SIZE];
  for (int i = 0; i < SIZE; i ++) {
    image[i] = TRUEVELOCITY;
  }

  VelocityField<double> velocity(image, WIDTH, HEIGHT);

  double mean_traveltimeerr = 0.0;
  double mean_reltraveltimeerr = 0.0;
  double mean_velerr = 0.0;
  double mean_relvelerr = 0.0;
  int meann = 0;
  double delta;

  for (size_t i = 0; i < STATIONS.size() - 1; i ++) {
    
    TravelTimeField<CartesianKm, double> tt(CartesianKm(LONMIN, LATMIN),
 					    CartesianKm(LONMAX, LATMAX),
					    CARTESIAN_STATIONS[i],
					    velocity,
					    1);
    
    tt.construct_traveltime_field();

    for (size_t j = i + 1; j < STATIONS.size(); j ++) {

      double traveltime = tt.get_traveltime_cubic(CARTESIAN_STATIONS[j]);
      double distkm = CartesianKm::distance_km(CARTESIAN_STATIONS[i], CARTESIAN_STATIONS[j]);
      double expected_traveltime = distkm/TRUEVELOCITY;
      double vel = distkm/traveltime;
      
      double traveltimeerr = fabs(traveltime - expected_traveltime);
      double reltraveltimeerr = traveltimeerr/expected_traveltime;

      double velerr = fabs(vel - TRUEVELOCITY);
      double relvelerr = velerr/TRUEVELOCITY;

      // printf("%5s %5s : %10.6fkm : %10.6fs %10.6fs : %10.6fkm/s %10.6fkm/s\n",
      //  	     STATIONNAMES[i].c_str(), STATIONNAMES[j].c_str(),
      //  	     distkm,
      //  	     traveltime, distkm/TRUEVELOCITY,
      //  	     distkm/traveltime, TRUEVELOCITY);

      meann ++;
      delta = traveltimeerr - mean_traveltimeerr;
      mean_traveltimeerr += delta/(double)meann;
      
      delta = reltraveltimeerr - mean_reltraveltimeerr;
      mean_reltraveltimeerr += delta/(double)meann;

      delta = velerr - mean_velerr;
      mean_velerr += delta/(double)meann;
      
      delta = relvelerr - mean_relvelerr;
      mean_relvelerr += delta/(double)meann;
    }
  }

  printf("%3d Traveltime %10.6f %10.6f Velocity %10.6f %10.6f\n",
	 meann,
	 mean_traveltimeerr, mean_reltraveltimeerr,
	 mean_velerr, mean_relvelerr);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
      

						     

    
  
  



  
