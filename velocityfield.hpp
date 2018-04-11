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
#ifndef velocityfield_hpp
#define velocityfield_hpp

template
<
  typename real
>
class VelocityField {
public:

  VelocityField(size_t _width,
		size_t _height) :
    field(nullptr),
    width(_width),
    height(_height),
    dx(1.0/(double)width),
    dy(1.0/(double)height)
  {
  }

  VelocityField(const real *_field,
		size_t _width,
		size_t _height) :
    field(_field),
    width(_width),
    height(_height),
    dx(1.0/(double)width),
    dy(1.0/(double)height)
  {
  }

  void set_field(const real *_field)
  {
    field = _field;
  }

  real lerp(double nx, double ny) const
  {
    //
    // Bilinear interpolation using normalized coordinates
    //
    int ix;
    double alphax;

    clamped_interpolation_coordinates(nx, dx, width, ix, alphax);

    int iy;
    double alphay;

    clamped_interpolation_coordinates(ny, dy, height, iy, alphay);

    real top = (1.0 - alphax)*field[iy * width + ix] + alphax * field[iy * width + ix + 1];
    real bot = (1.0 - alphax)*field[(iy + 1) * width + ix] + alphax * field[(iy + 1) * width + ix + 1];

    return (1.0 - alphay)*top + alphay * bot;
  }

  void clamped_interpolation_coordinates(double nx, double dx, size_t sizex, int &ix, double &alphax) const
  {
    ix = (int)((double)sizex * (nx - dx/2.0));
    if (ix < 0) {
      ix = 0;
      alphax = 0.0;
    } else if (ix >= (int)(sizex - 1)) {
      ix = sizex - 2;
      alphax = 1.0;
    } else {
      alphax = (nx - ((double)ix + 0.5) * dx)/dx;
      if (alphax < 0.0) {
	alphax = 0.0;
      }
    }
  }

  size_t get_width() const
  {
    return width;
  }

  size_t get_height() const
  {
    return height;
  }

  void dump(FILE *fp) const
  {
    for (size_t j = 0; j < height; j ++) {
      for (size_t i = 0; i < width; i ++) {
	fprintf(fp, "%7.3f ", field[j * width + i]);
      }
      fprintf(fp, "\n");
    }
  }

  const real *field;
  
private:

  size_t width;
  size_t height;

  double dx;
  double dy;
  
};


#endif // velocityfield_hpp

  
