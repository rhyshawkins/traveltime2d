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
#include "velocityweights.hpp"

template
<
  typename coordinate,
  typename real
> class TravelTimeField;

template
<
  typename coordinate,
  typename real
>
class TravelTimeNode {
public:

  typedef VelocityField<real> velocity_t;

  struct TravelTimeWeight {
    int n;
    TravelTimeNode *neighbors[4];
    real weight[4];
    real vl_weight;

    TravelTimeNode *linked_node;
  };

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

    weight.n = 0;
    weight.vl_weight = 0.0;
    weight.linked_node = nullptr;

    vweights.reset();

    vweights.n = 4;
    vf.lerp_weights(nx, ny,
		    vweights.indices[0], vweights.weights[0],
		    vweights.indices[1], vweights.weights[1],
		    vweights.indices[2], vweights.weights[2],
		    vweights.indices[3], vweights.weights[3]);

    
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

  void initializeT(double _T,
		   TravelTimeNode *_node)
  {
    //
    // Initialize a fixed velocity node
    //
    T = _T;
    weight.linked_node = _node;

    state = STATE_FIXED;
  }

  void initializeT(const velocity_t &vf, double distkm)
  {
    //
    // Initialize a fixed velocity node
    //
    T = distkm/V;

    weight.n = 0;
    weight.vl_weight = -distkm/(V*V);
    weight.linked_node = nullptr;

    // printf("Initial weight: %10.6f\n", weight.vl_weight);
    
    state = STATE_FIXED;
  }

  real computeTdirectionderivative(TravelTimeNode *h1,
				   TravelTimeNode *h2,
				   const real &dist,
				   TravelTimeWeight &weight)
  {
    if (h2 == nullptr || firstorder) {
      weight.n = 1;
      weight.neighbors[0] = h1;
      weight.weight[0] = 1.0;
      weight.vl_weight = -dist/(V*V);
      
      return h1->T + dist/V;
      
    } else {
      real T = ((4.0*h1->T - h2->T) + 2.0*dist/V)/3.0;

      if (T < h1->T) {
	weight.n = 1;
	weight.neighbors[0] = h1;
	weight.weight[0] = 1.0;
	weight.vl_weight = -dist/(V*V);

	T = h1->T + dist/V;

      } else {
	weight.n = 2;

	weight.neighbors[0] = h1;
	weight.weight[0] = 4.0/3.0;

	weight.neighbors[1] = h2;
	weight.weight[1] = -1.0/3.0;
	
	weight.vl_weight = -2.0*dist/(3.0*V*V);
      }
      
      return T;
    }
  }

  real computeTdualderivative(TravelTimeNode *h1, TravelTimeNode *h2, const real &hdist,
			      TravelTimeNode *v1, TravelTimeNode *v2, const real &vdist,
			      TravelTimeWeight &weight)
  {
    real Tx, Dx;
    real Ty, Dy;
    bool dualx, dualy;

    if (firstorder || h2 == nullptr) {
      Tx = h1->T;
      Dx = hdist;
      dualx = false;
    } else {
      Tx = (4.0*h1->T - h2->T)/3.0;
      Dx = 2.0*hdist/3.0;
      dualx = true;
    }

    if (firstorder || v2 == nullptr) {
      Ty = v1->T;
      Dy = vdist;
      dualy = false;
    } else {
      Ty = (4.0*v1->T - v2->T)/3.0;
      Dy = 2.0*vdist/3.0;
      dualy = true;
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
	TravelTimeWeight weightA, weightB;
	
	Ta = computeTdirectionderivative(h1,
					 h2,
					 hdist,
					 weightA);
	Tb = computeTdirectionderivative(v1,
					 v2,
					 vdist,
					 weightB);
	
	if (Ta < Tb) {
	  weight = weightA;
	  return Ta;
	} else {
	  weight = weightB;
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
	TravelTimeWeight weightA, weightB;

	Ta = computeTdirectionderivative(h1,
					 h2,
					 hdist,
					 weightA);
	
	Tb = computeTdirectionderivative(v1,
					 v2,
					 vdist,
					 weightB);
	
	if (Ta < Tb) {
	  weight = weightA;
	  return Ta;
	} else {
	  weight = weightB;
	  return Tb;
	}
	
	// throw TRAVELTIMEEXCEPTION("Causality violated %f %f (%f %f %g %g %g %g) (%f %f) %f", Ta, Tb,
	// 			    h1->T, v1->T,
	// 			    Ta - h1->T, Ta - v1->T,
	// 			    Tb - h1->T, Tb - v1->T,
	// 			    hdist, vdist,
	// 			    V);
	
      }

      real denominator = 2.0*(Tb * (1.0/Dx2 + 1.0/Dy2)) - 2.0*((Tx/Dx2 + Ty/Dy2));
      real wh = (2.0*(Tb - Tx)/Dx2)/denominator;
      real wv = (2.0*(Tb - Ty)/Dy2)/denominator;

      if (dualx) {
	weight.neighbors[0] = h1;
	weight.neighbors[1] = h2;

	weight.weight[0] = wh * 4.0/3.0;
	weight.weight[1] = -wh/3.0;

	if (dualy) {

	  weight.neighbors[2] = v1;
	  weight.neighbors[3] = v2;
	  
	  weight.weight[2] = wv * 4.0/3.0;
	  weight.weight[3] = -wv/3.0;

	  weight.n = 4;
	  
	} else {
	  weight.neighbors[2] = v1;
	  weight.weight[2] = wv;
	  weight.n = 3;
	}
	
      } else {
	
	weight.neighbors[0] = h1;
	weight.weight[0] = wh;

	if (dualy) {

	  weight.neighbors[1] = v1;
	  weight.neighbors[2] = v2;
	  
	  weight.weight[1] = wv * 4.0/3.0;
	  weight.weight[2] = -wv/3.0;

	  weight.n = 3;
	  
	} else {
	  
	  weight.neighbors[1] = v1;
	  weight.weight[1] = wv;
	  weight.n = 2;
	  
	}
	
      }

      weight.vl_weight = (-2.0/(V*V*V))/denominator;
      
      return Tb;
      
    } else {
      
      // Choose smallest (Ta always less than Tb)
      real denominator = 2.0*(Ta * (1.0/Dx2 + 1.0/Dy2)) - 2.0*((Tx/Dx2 + Ty/Dy2));
      real wh = (2.0*(Ta - Tx)/Dx2)/denominator;
      real wv = (2.0*(Ta - Ty)/Dy2)/denominator;

      if (dualx) {
	weight.neighbors[0] = h1;
	weight.neighbors[1] = h2;

	weight.weight[0] = wh * 4.0/3.0;
	weight.weight[1] = -wh/3.0;

	if (dualy) {

	  weight.neighbors[2] = v1;
	  weight.neighbors[3] = v2;
	  
	  weight.weight[2] = wv * 4.0/3.0;
	  weight.weight[3] = -wv/3.0;

	  weight.n = 4;
	  
	} else {
	  weight.neighbors[2] = v1;
	  weight.weight[2] = wv;
	  weight.n = 3;
	}
	
      } else {
	
	weight.neighbors[0] = h1;
	weight.weight[0] = wh;

	if (dualy) {

	  weight.neighbors[1] = v1;
	  weight.neighbors[2] = v2;
	  
	  weight.weight[1] = wv * 4.0/3.0;
	  weight.weight[2] = -wv/3.0;

	  weight.n = 3;
	  
	} else {
	  
	  weight.neighbors[1] = v1;
	  weight.weight[1] = wv;
	  weight.n = 2;
	  
	}
	
      }

      weight.vl_weight = (-2.0/(V*V*V))/denominator;
      
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
				    vertical, vertical_second, vdist,
				    weight);

    } else if (horizontal) {

      return computeTdirectionderivative(horizontal,
					 horizontal_second,
					 hdist,
					 weight);

    } else if (vertical) {

      return computeTdirectionderivative(vertical,
					 vertical_second,
					 vdist,
					 weight);

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

  void back_project(TravelTimeField<coordinate, real> *parent)
  {
    if (vweights.dirty) {

      if (weight.linked_node == nullptr) {
	
	for (int i = 0; i < 4; i ++) {
	  vweights.weights[i] *= weight.vl_weight;
	  // if (vweights.indices[i] == 8255) {
	  //   printf("Setting: %5d %10.6f\n", vweights.indices[i], vweights.weights[i]);
	  // }
	}
	
	for (int i = 0; i < weight.n; i ++) {
	  
	  weight.neighbors[i]->back_project(parent);
	  
	  vweights.merge(weight.neighbors[i]->vweights,
			 weight.weight[i]);
	  
	}

	if (!vweights.validate()) {
	  throw TRAVELTIMEEXCEPTION("Invalid vweights");
	}

      } else {

	weight.linked_node->back_project(parent->subtraveltime);

	//
	// Remap weights from sub field to proper indices into
	// image
	//

	vweights.n = 0;
	
	for (int i = 0; i < weight.linked_node->vweights.n; i ++) {

	  auto it = parent->subfield->interpolators[weight.linked_node->vweights.indices[i]];

	  it->backproject_remap(vweights,
				weight.linked_node->vweights.weights[i]);
	}
      }
      
      vweights.dirty = false;
    }
  }
  
  double nx, ny;
  state_t state;
  int order;
  bool dirty;
  real T;
  real V;

  TravelTimeWeight weight;
  VelocityWeights<real> vweights;
  
  Neighbor neighbors[4];

  bool firstorder;
};

#endif // traveltimenode_hpp

