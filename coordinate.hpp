//
//    TravelTime2d : A library for computing travel times for 2D surface
//    wave propagation problems with the Fast Marching method.
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once
#ifndef coordinate_hpp
#define coordinate_hpp

#include <math.h>

template
<
  //
  // I'm assuming we will always do 2D on surface of a sphere so radius is constant.
  // This is just a way to encode the radius as a constant and override at compile time.
  //
  int radius_numerator = 6371,
  int radius_denominator = 1
>
class LonLat {
public:

  LonLat() :
    lon(0.0),
    lat(0.0)
  {
  }
  
  LonLat(double _lon, double _lat) :
    lon(_lon),
    lat(_lat)
  {
  }

  double X() const {
    return lon;
  }

  double Y() const {
    return lat;
  }

  static LonLat mesh_coordinate(const LonLat &cmin,
				const LonLat &cmax,
				size_t width, // Velocity field dimensions
				size_t height,
				int i,
				int j)
  {
    double dlon = (cmax.lon - cmin.lon)/(double)(width);
    double lon = cmin.lon - dlon/2.0 + dlon * i;

    double dlat = (cmax.lat - cmin.lat)/(double)(height);
    double lat = cmin.lat - dlat/2.0 + dlat * j;

    return LonLat(lon, lat);
  }

  static LonLat midpoint(const LonLat &a, const LonLat &b)
  {
    return LonLat((a.lon + b.lon)/2.0,
		  (a.lat + b.lat)/2.0);
  }

  static void normalized_coordinate(const LonLat &cmin,
				    const LonLat &cmax,
				    const LonLat &c,
				    double &nx,
				    double &ny)
  {
    nx = (c.lon - cmin.lon)/(cmax.lon - cmin.lon);
    ny = (c.lat - cmin.lat)/(cmax.lat - cmin.lat);
  }

  static LonLat from_normalized_coordinate(const LonLat &cmin,
					   const LonLat &cmax,
					   double nx, double ny)
  {
    return LonLat(cmin.lon + nx * (cmax.lon - cmin.lon),
		  cmin.lat + ny * (cmax.lat - cmin.lat));
  }

  static double distance_km(const LonLat &a,
			    const LonLat &b)
  {
    return distance_radians(a, b) * (double)radius_numerator/(double)radius_denominator;
  }

  static double distance_radians(const LonLat &a,
				 const LonLat &b)
  {
    double lon1 = M_PI/180.0 * a.lon;
    double lat1 = M_PI/180.0 * a.lat;

    double lon2 = M_PI/180.0 * b.lon;
    double lat2 = M_PI/180.0 * b.lat;

    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;

    double hsinlat = sin(dlat/2.0);
    double hsinlon = sin(dlon/2.0);
    
    return 2.0 * asin(sqrt(hsinlat * hsinlat + cos(lat1)*cos(lat2)*(hsinlon*hsinlon)));
  }

  double lon, lat;
};

class CartesianKm {
public:

  CartesianKm() :
    x(0.0),
    y(0.0)
  {
  }
  
  CartesianKm(double _x, double _y) :
    x(_x),
    y(_y)
  {
  }

  double X() const {
    return x;
  }

  double Y() const {
    return y;
  }

  static CartesianKm mesh_coordinate(const CartesianKm &cmin,
				     const CartesianKm &cmax,
				     size_t vfwidth, 
				     size_t vfheight,
				     int i,
				     int j)
  {
    double dx = (cmax.x - cmin.x)/(double)(vfwidth);
    double x = cmin.x - dx/2.0 + dx * i;

    double dy = (cmax.y - cmin.y)/(double)(vfheight);
    double y = cmin.y - dy/2.0 + dy * j;

    return CartesianKm(x, y);
  }

  static CartesianKm midpoint(const CartesianKm &a, const CartesianKm &b)
  {
    return CartesianKm((a.x + b.x)/2.0,
		       (a.y + b.y)/2.0);
  }

  static void normalized_coordinate(const CartesianKm &cmin,
				    const CartesianKm &cmax,
				    const CartesianKm &c,
				    double &nx,
				    double &ny)
  {
    nx = (c.x - cmin.x)/(cmax.x - cmin.x);
    ny = (c.y - cmin.y)/(cmax.y - cmin.y);
  }
  
  static CartesianKm from_normalized_coordinate(const CartesianKm &cmin,
						const CartesianKm &cmax,
						double nx, double ny)
  {
    return CartesianKm(cmin.x + nx * (cmax.x - cmin.x),
		       cmin.y + ny * (cmax.y - cmin.y));
  }


  static double distance_km(const CartesianKm &a,
			    const CartesianKm &b)
  {
    double dx = b.x - a.x;
    double dy = b.y - a.y;

    return sqrt(dx*dx + dy*dy);
  }

  double x, y;
};

#endif // coordinate.hpp
