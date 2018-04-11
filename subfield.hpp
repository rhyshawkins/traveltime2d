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
#ifndef subfield_hpp
#define subfield_hpp

#include "traveltimeexception.hpp"

template
<
  typename real
>
class Subfield {
public:

  Subfield(int source_width,
	   int source_height,
	   int source_i,
	   int source_j,
	   int _width,
	   int _height) :
    width(_width),
    height(_height),
    size(_width * _height),
    field(new real[_width * _height]),
    interpolators(new interpolant*[_width * _height])
  {
    // printf("Subfield: %d %d: %d %d %d %d\n",
    // 	   source_width,
    // 	   source_height,
    // 	   source_i,
    // 	   source_j,
    // 	   _width,
    // 	   _height);
    
    int source_i_end = source_i + (width + 1)/2;
    if (source_i_end >= source_width) {
      //throw TRAVELTIMEEXCEPTION("Source overflow horizontally %d .. %d (%d)", source_i, source_i_end, source_width);
      source_i_end = source_width - 1;
    }

    int source_j_end = source_j + (height + 1)/2;
    if (source_j_end >= source_height) {
      //throw TRAVELTIMEEXCEPTION("Source overflow vertically");
      source_j_end = source_height - 1;
    }

    for (int j = 0; j < height; j ++) {

      for (int i = 0; i < width; i ++) {

	if (j % 2 == 0) {

	  if (i % 2 == 0) {

	    interpolators[j * width + i] = new single_interpolant((source_j + j/2)*source_width + source_i + i/2);

	  } else {

	    interpolators[j * width + i] = new linear_interpolant((source_j + j/2)*source_width + source_i + i/2,
								  (source_j + j/2)*source_width + source_i + i/2 + 1);
	  }

	} else {

	  if (i % 2 == 0) {
	    interpolators[j * width + i] = new linear_interpolant((source_j + j/2)*source_width + source_i + i/2,
								  (source_j + j/2 + 1)*source_width + source_i + i/2);
	  } else {
	    interpolators[j * width + i] = new bilinear_interpolant((source_j + j/2)*source_width + source_i + i/2,
								    (source_j + j/2)*source_width + source_i + i/2 + 1,
								    (source_j + j/2 + 1)*source_width + source_i + i/2,
								    (source_j + j/2 + 1)*source_width + source_i + i/2 + 1);

	  }

	}
      }
    }
  }

  void update(const real *source)
  {
    for (int i = 0; i < size; i ++) {
      field[i] = interpolators[i]->interpolate(source);
      if (field[i] <= 0.0) {
	interpolators[i]->dump_source();
	throw TRAVELTIMEEXCEPTION("Invalid field value %f (%d/%d)\n", field[i], i, size);
      }
    }
  }

  class interpolant {
  public:
    interpolant()
    {
    }
    virtual ~interpolant()
    {
    }

    virtual real interpolate(const real *source) const = 0;

    virtual void dump_source() const = 0;
  };

  class single_interpolant : public interpolant {
  public:
    single_interpolant(int _i) :
      i(_i)
    {
      if (i < 0) {
	throw TRAVELTIMEEXCEPTION("Invalid index: %d", _i);
      }
    }
    
    virtual real interpolate(const real *source) const
    {
      return source[i];
    }

    virtual void dump_source() const
    {
      printf("Single: %d\n", i);
    }
    
    int i;
  };

  class linear_interpolant : public interpolant {
  public:
    linear_interpolant(int _i, int _j) :
      i(_i),
      j(_j)
    {
    }
    
    virtual real interpolate(const real *source) const
    {
      return (source[i] + source[j])/2.0;
    }

    virtual void dump_source() const
    {
      printf("Linear: %d %d\n", i, j);
    }

    int i;
    int j;
  };

  class bilinear_interpolant : public interpolant {
  public:
    bilinear_interpolant(int _i, int _j, int _k, int _l) :
      i(_i),
      j(_j),
      k(_k),
      l(_l)
    {
    }
    
    virtual real interpolate(const real *source) const
    {
      return (source[i] + source[j] + source[k] + source[l])/4.0;
    }

    virtual void dump_source() const
    {
      printf("Bilinear: %d %d %d %d\n", i, j, k, l);
    }
    int i;
    int j;
    int k;
    int l;
    
  };

  int width;
  int height;
  int size;
  real *field;
  interpolant **interpolators;
};

#endif // subfield_hpp
