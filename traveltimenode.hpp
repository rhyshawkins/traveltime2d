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
#ifndef traveltimenode_hpp
#define traveltimenode_hpp

#include <stdio.h>

#include "velocityfield.hpp"
#include "traveltimeexception.hpp"

template
<
  typename real
>
class TravelTimeNode {
public:

  typedef VelocityField<real> velocity_t;
  
  static constexpr real INVALID_T = -1.0;
  
  struct Neighbor {

    Neighbor() :
      neighbor(nullptr),
      distkm(0.0)
    {
    }

    Neighbor(TravelTimeNode &_neighbor,
	     double _distkm) :
      neighbor(&_neighbor),
      distkm(_distkm)
    {
    }
    
    TravelTimeNode *neighbor;
    double distkm;
  };

  typedef enum {
    STATE_FIXED,
    STATE_NEAR,
    STATE_FAR
  } state_t;

  typedef enum {
    NEIGHBOR_LEFT,
    NEIGHBOR_RIGHT,
    NEIGHBOR_TOP,
    NEIGHBOR_BOTTOM
  } neighbor_direction_t;
    
  TravelTimeNode() :
    nx(0.0),
    ny(0.0),
    state(STATE_FAR),
    order(-1),
    dirty(true),
    T(INVALID_T),
    V(0.0),
    firstorder(false)
  {
  }

  void reset(const velocity_t &vf)
  {
    state = STATE_FAR;
    T = INVALID_T;
    order = -1;
    V = vf.lerp(nx, ny);
    
    if (V <= 0.0) {
      vf.dump(stderr);
      throw TRAVELTIMEEXCEPTION("Invalid velocity: %f", V);
    }
  }

  void set_order(int o)
  {
    order = o;
  }

  int get_order() const
  {
    return order;
  }
  
  void set_normalized_coordinates(double _nx, double _ny)
  {
    nx = _nx;
    ny = _ny;
  }

  void add_top_neighbor(TravelTimeNode &other,
			double distkm)
  {
    neighbors[NEIGHBOR_TOP] = Neighbor(other, distkm);
    other.neighbors[NEIGHBOR_BOTTOM] = Neighbor(*this, distkm);
  }
  
  void add_left_neighbor(TravelTimeNode &other,
			 double distkm)
  {
    neighbors[NEIGHBOR_LEFT] = Neighbor(other, distkm);
    other.neighbors[NEIGHBOR_RIGHT] = Neighbor(*this, distkm);
  }

  void initializeT(double _T)
  {
    //
    // Initialize a fixed velocity node
    //
    T = _T;
    state = STATE_FIXED;
  }

  void initializeT(const velocity_t &vf, double distkm)
  {
    //
    // Initialize a fixed velocity node
    //
    T = distkm/V;
    state = STATE_FIXED;
  }

  real computeTdirectionderivative(TravelTimeNode *h1,
				   TravelTimeNode *h2,
				   const real &dist)
  {
    if (h2 == nullptr || firstorder) {
      return h1->T + dist/V;
    } else {
      real T = ((4.0*h1->T - h2->T) + 2.0*dist/V)/3.0;
      if (T < h1->T) {
	T = h1->T + dist/V;
      }
      return T;
    }
  }

  real computeTdualderivative(TravelTimeNode *h1, TravelTimeNode *h2, const real &hdist,
			      TravelTimeNode *v1, TravelTimeNode *v2, const real &vdist)
  {
    real Tx, Dx;
    real Ty, Dy;

    if (firstorder || h2 == nullptr) {
      Tx = h1->T;
      Dx = hdist;
    } else {
      Tx = (4.0*h1->T - h2->T)/3.0;
      Dx = 2.0*hdist/3.0;
    }

    if (firstorder || v2 == nullptr) {
      Ty = v1->T;
      Dy = vdist;
    } else {
      Ty = (4.0*v1->T - v2->T)/3.0;
      Dy = 2.0*vdist/3.0;
    }

    real Dx2 = Dx*Dx;
    real Dy2 = Dy*Dy;

    real Ta, Tb;

    solve_quadratic(1.0/Dx2 + 1.0/Dy2,
		    -2.0*(Tx/Dx2 + Ty/Dy2),
		    Tx*Tx/Dx2 + Ty*Ty/Dy2 - 1.0/(V*V),
		    Ta, Tb);
    
    if (Tb < 0.0) {
      if (Ta < 0.0) {
	Ta = computeTdirectionderivative(h1,
					 h2,
					 hdist);
	Tb = computeTdirectionderivative(v1,
					 v2,
					 vdist);
	
	if (Ta < Tb) {
	  return Ta;
	} else {
	  return Tb;
	}
	
      } else {
	if (Ta < h1->T && Ta < v1->T) {
	  throw TRAVELTIMEEXCEPTION("Causality violated %f %f (%f %f) (%f %f) %f", Ta, Tb,
				    h1->T, v1->T,
				    hdist, vdist,
				    V);
	  
	}
	
	return Ta;
      }
    }

    //
    // Two solutions Ta will be less that Tb always
    //
    if (Ta < h1->T || Ta < v1->T) {
      
      if (Tb < h1->T || Tb < v1->T) {
	
	//
	// Numerical precision can mean that causaility fails for uniform velocity fields
	// so revert to best orthogonal travel time in this case.
	//
	Ta = computeTdirectionderivative(h1,
					 h2,
					 hdist);
	
	Tb = computeTdirectionderivative(v1,
					 v2,
					 vdist);
	
	if (Ta < Tb) {
	  return Ta;
	} else {
	  return Tb;
	}
	
	// throw TRAVELTIMEEXCEPTION("Causality violated %f %f (%f %f %g %g %g %g) (%f %f) %f", Ta, Tb,
	// 			    h1->T, v1->T,
	// 			    Ta - h1->T, Ta - v1->T,
	// 			    Tb - h1->T, Tb - v1->T,
	// 			    hdist, vdist,
	// 			    V);
	
      }
      
      return Tb;
      
    } else {
      
      // Choose smallest (Ta always less than Tb)
      return Ta;
      
    }

    throw TRAVELTIMEEXCEPTION("Unreachable");
    return INVALID_T;
  }

  
  real computeT(const velocity_t &vf)
  {
    //
    // First determine which horizontal directions are available
    //
    TravelTimeNode *horizontal = nullptr;
    TravelTimeNode *horizontal_second = nullptr;
    double hdist = 0.0;

    if (neighbors[NEIGHBOR_LEFT].neighbor != nullptr &&
	neighbors[NEIGHBOR_LEFT].neighbor->state == STATE_FIXED) {

      horizontal = neighbors[NEIGHBOR_LEFT].neighbor;
      hdist = neighbors[NEIGHBOR_LEFT].distkm;

      if (horizontal->neighbors[NEIGHBOR_LEFT].neighbor != nullptr &&
	  horizontal->neighbors[NEIGHBOR_LEFT].neighbor->state == STATE_FIXED) {
	horizontal_second = horizontal->neighbors[NEIGHBOR_LEFT].neighbor;
      }
    }

    if (neighbors[NEIGHBOR_RIGHT].neighbor != nullptr &&
	neighbors[NEIGHBOR_RIGHT].neighbor->state == STATE_FIXED) {

      if (horizontal == nullptr ||
	  (neighbors[NEIGHBOR_RIGHT].neighbor->T <
	   neighbors[NEIGHBOR_LEFT].neighbor->T)) {
	
	horizontal = neighbors[NEIGHBOR_RIGHT].neighbor;
	hdist = neighbors[NEIGHBOR_RIGHT].distkm;
	
	if (horizontal->neighbors[NEIGHBOR_RIGHT].neighbor != nullptr &&
	    horizontal->neighbors[NEIGHBOR_RIGHT].neighbor->state == STATE_FIXED) {
	  horizontal_second = horizontal->neighbors[NEIGHBOR_RIGHT].neighbor;
	} else {
	  // Incase the second horizontal was set for the left direction
	  horizontal_second = nullptr;
	}
	
      }   
    }
	
    //
    // Second determine which vertical directions are available
    //
    TravelTimeNode *vertical = nullptr;
    TravelTimeNode *vertical_second = nullptr;
    double vdist = 0.0;
    
    if (neighbors[NEIGHBOR_TOP].neighbor != nullptr &&
	neighbors[NEIGHBOR_TOP].neighbor->state == STATE_FIXED) {

      vertical = neighbors[NEIGHBOR_TOP].neighbor;
      vdist = neighbors[NEIGHBOR_TOP].distkm;

      if (vertical->neighbors[NEIGHBOR_TOP].neighbor != nullptr &&
	  vertical->neighbors[NEIGHBOR_TOP].neighbor->state == STATE_FIXED) {
	vertical_second = vertical->neighbors[NEIGHBOR_TOP].neighbor;
      }
      
    }

    if (neighbors[NEIGHBOR_BOTTOM].neighbor != nullptr &&
	neighbors[NEIGHBOR_BOTTOM].neighbor->state == STATE_FIXED) {

      if (vertical == nullptr ||
	  (neighbors[NEIGHBOR_BOTTOM].neighbor->T <
	   neighbors[NEIGHBOR_TOP].neighbor->T)) {
	
	vertical = neighbors[NEIGHBOR_BOTTOM].neighbor;
	vdist = neighbors[NEIGHBOR_BOTTOM].distkm;

	if (vertical->neighbors[NEIGHBOR_BOTTOM].neighbor != nullptr &&
	    vertical->neighbors[NEIGHBOR_BOTTOM].neighbor->state == STATE_FIXED) {
	  vertical_second = vertical->neighbors[NEIGHBOR_BOTTOM].neighbor;
	} else {
	  // Incase the second vertical was set for the left direction
	  vertical_second = nullptr;
	}
	
      }   
    }

    if (horizontal && vertical) {

      return computeTdualderivative(horizontal, horizontal_second, hdist,
				    vertical, vertical_second, vdist);

      real Ta, Tb;
      //
      // Solve quadratic
      //

      
      solve_quadratic(vdist*vdist + hdist*hdist,
		      -2.0 * (vdist*vdist * horizontal->T + hdist*hdist * vertical->T),
		      vdist*vdist * horizontal->T*horizontal->T +
		      hdist*hdist * vertical->T*vertical->T -
		      (hdist*hdist*vdist*vdist)/(V*V),
		      Ta, Tb);

      // printf("%f %f %f %f: %f %f\n", horizontal->T, hdist, vertical->T, vdist, Ta, Tb);
      
      if (Tb < 0.0) {
	if (Ta < 0.0) {
	  Ta = computeTdirectionderivative(horizontal,
					   horizontal_second,
					   hdist);
	  Tb = computeTdirectionderivative(vertical,
					   vertical_second,
					   vdist);

	  if (Ta < Tb) {
	    return Ta;
	  } else {
	    return Tb;
	  }
	} else {
	  if (Ta < horizontal->T && Ta < vertical->T) {
	    throw TRAVELTIMEEXCEPTION("Causality violated %f %f (%f %f) (%f %f) %f", Ta, Tb,
				      horizontal->T, vertical->T,
				      hdist, vdist,
				      V);
	    
	  }
	  
	  return Ta;
	}
      }

      //
      // Two solutions Ta will be less that Tb always
      //
      if (Ta < horizontal->T || Ta < vertical->T) {
	
	if (Tb < horizontal->T || Tb < vertical->T) {

	  //
	  // Numerical precision can mean that causaility fails for uniform velocity fields
	  // so revert to best orthogonal travel time in this case.
	  //
	  Ta = computeTdirectionderivative(horizontal,
					   horizontal_second,
					   hdist);

	  Tb = computeTdirectionderivative(vertical,
					   vertical_second,
					   vdist);

	  if (Ta < Tb) {
	    return Ta;
	  } else {
	    return Tb;
	  }

	  // throw TRAVELTIMEEXCEPTION("Causality violated %f %f (%f %f %g %g %g %g) (%f %f) %f", Ta, Tb,
	  // 			    horizontal->T, vertical->T,
	  // 			    Ta - horizontal->T, Ta - vertical->T,
	  // 			    Tb - horizontal->T, Tb - vertical->T,
	  // 			    hdist, vdist,
	  // 			    V);

	}

	return Tb;

      } else {

	// Choose smallest (Ta always less than Tb)
	return Ta;

      }

    } else if (horizontal) {

      return computeTdirectionderivative(horizontal,
					 horizontal_second,
					 hdist);

    } else if (vertical) {

      return computeTdirectionderivative(vertical,
					 vertical_second,
					 vdist);

    } else {

      //
      // Error
      //
      printf("No horizontal or vertical fixed nodes\n");
      return INVALID_T;

    }
  }

  void solve_quadratic(real a, real b, real c, real &Ta, real &Tb)
  {
    real d = b*b - 4.0 * a * c;
    
    if (d < 0.0) {
      Ta = INVALID_T;
      Tb = INVALID_T;
    } else if (d == 0.0) {
      Ta = -b/(2.0 * a);
      Tb = INVALID_T;
    } else {
      Ta = (-b - sqrt(d))/(2.0 * a);
      Tb = (-b + sqrt(d))/(2.0 * a);
    }      
  }

  double nx, ny;
  state_t state;
  int order;
  bool dirty;
  real T;
  real V;
  
  Neighbor neighbors[4];

  bool firstorder;
};

#endif // traveltimenode_hpp

