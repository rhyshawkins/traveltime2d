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
#ifndef velocityfieldref_hpp
#define velocityfieldref_hpp

template
<
  typename real
>
class VelocityFieldRef {
public:

  VelocityFieldRef() :
    dirty(true),
    value(0.0)
  {
    for (int i = 0; i < 4; i ++) {
      ref[0] = nullptr;
      scale[0] = nullptr;
    }
  }

  real evaluate()
  {
    if (dirty) {
      value =
	(scale[0] * (*ref[0])) + 
	(scale[1] * (*ref[1])) + 
	(scale[2] * (*ref[2])) + 
	(scale[3] * (*ref[3]));

      dirty = false;
    }

    return value;
  }

  real *ref[4];
  real scale[4];

  bool dirty;
  real value;
};

#endif // velocityfieldref_hpp

