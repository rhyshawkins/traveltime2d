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
#ifndef traveltimefield_hpp
#define traveltimefield_hpp

#include <set>
#include <vector>

#include "coordinate.hpp"
#include "velocityfield.hpp"
#include "traveltimenode.hpp"
#include "circularbuffer.hpp"
#include "traveltimeexception.hpp"
#include "subfield.hpp"

template
<
  typename coordinate,
  typename real
>
class TravelTimeField {
public:

  typedef TravelTimeNode<real> node_t;
  typedef VelocityField<real> velocity_t;

  static constexpr real LARGE_T = 1.0e12;
  static constexpr int DIRTY_SIZE = 32;

  TravelTimeField(const coordinate &_cmin,
		  const coordinate &_cmax,
		  const coordinate &_hypo,
		  velocity_t &_velocity,
		  int _refinement) :
    cmin(_cmin),
    cmax(_cmax),
    hypo(_hypo),
    velocity(_velocity),
    width(velocity.get_width() + 2),
    height(velocity.get_height() + 2),
    refinement(_refinement),
    near(2*(velocity.get_width() + velocity.get_height())),
    dirty_count(0),
    nodes(nullptr),
    dx(1.0/(real)velocity.get_width()),
    dy(1.0/(real)velocity.get_height()),
    nminx(-0.5*dx),
    nmaxx(1.0 + 0.5*dx),
    nminy(-0.5*dy),
    nmaxy(1.0 + 0.5*dy)
  {
    if (cmin.X() >= cmax.X() ||
	cmin.Y() >= cmax.Y()) {
      throw TRAVELTIMEEXCEPTION("Invalid boundary coordinates");
    }

    
    nodes = new node_t[width * height];

    //
    // Initialize the coordinates of the mesh nodes
    //
    for (size_t j = 0; j < height; j ++) {

      for (size_t i = 0; i < width; i ++) {

	double nx, ny;
	
	coordinate::normalized_coordinate(cmin, cmax,
					  coordinate::mesh_coordinate(cmin, cmax,
								      width - 2, height - 2,
								      i, j),
					  nx, ny);
	
	nodes[j * width + i].set_normalized_coordinates(nx, ny);

      }
    }

    //
    // Initialize neighbourhood information
    //
    for (size_t j = 0; j < height; j ++) {

      for (size_t i = 0; i < width; i ++) {

	int a = j * width + i;

	coordinate ca = coordinate::from_normalized_coordinate(cmin, cmax,
							       nodes[a].nx, nodes[a].ny);
	
	if (j > 0) {
	  int b = (j - 1) * width + i;	
	  coordinate cb = coordinate::from_normalized_coordinate(cmin, cmax,
								 nodes[b].nx, nodes[b].ny);
  
	  nodes[a].add_top_neighbor(nodes[b],
				    coordinate::distance_km(ca, cb));
	}

	if (i > 0) {

	  int b = j * width + i - 1;
	  coordinate cb = coordinate::from_normalized_coordinate(cmin, cmax,
								 nodes[b].nx, nodes[b].ny);
	  
	  nodes[a].add_left_neighbor(nodes[b],
				     coordinate::distance_km(ca, cb));
	}
      }
    }

    //
    // Construct Refinement grid
    //
    if (refinement > 0) {

      int bl, br, tl, tr;

      find_enclosing_nodes(hypo, bl, br, tl, tr);

      double dnx;
      double dny;

      dnx = nodes[br].nx - nodes[bl].nx;
      dny = nodes[tl].ny - nodes[bl].ny;
      
      double minnx, maxnx;

      minnx = nodes[bl].nx - 5.0/4.0 * dnx;
      maxnx = nodes[br].nx + 5.0/4.0 * dnx;

      double minny, maxny;
      
      minny = nodes[bl].ny - 5.0/4.0 * dny;
      maxny = nodes[tl].ny + 5.0/4.0 * dny;

      //
      // Create the sub image
      //
      int ni = bl % width;
      int nj = bl / width;

      if (ni < 2) {
	ni = 2;
      }
      if (nj < 2) {
	nj = 2;
      }

      //
      // Clamp ni/nj to valid ranges
      //

      // printf("Source i, j: %d %d (%d %d)\n", ni, nj, ni - 2, nj - 2);
      subfield = new Subfield<real>(width - 2,
				    height - 2,
				    ni - 2,
				    nj - 2,
				    7, 7);
      subvelocity = new velocity_t(subfield->field, 7, 7);

      coordinate subcmin = coordinate::from_normalized_coordinate(cmin, cmax, minnx, minny);
      coordinate subcmax = coordinate::from_normalized_coordinate(cmin, cmax, maxnx, maxny);

      //
      // Create the subfield
      //
      subtraveltime = new TravelTimeField(subcmin, subcmax, hypo, *subvelocity, refinement - 1);

      //
      // Create the node mappings
      //
      int o = bl - 1 - width;
      for (int j = 0; j < 4; j ++) {

	for (int i = 0; i < 4; i ++) {

	  subnodemappings.push_back(subnodemap_t(o + j * width + i, (2*j + 1) * 9 + 2*i + 1));
	  // printf("%d, %d : %d <- %d\n", i, j, (int)(o + j * width + i), (int)((2*j + 1) * 9 + 2*i + 1));
	}
      }

      
      //
      // Create the list of initial near nodes
      //
      int obl = bl - 1 - width;
      for (int i = 0; i < 4; i ++) {

	// Bottom row
	subnodenear.push_back(obl + i);

	// Top row
	subnodenear.push_back(obl + 3 * width + i);

      }

      for (int i = 1; i < 3; i ++) {

	// Left col
	subnodenear.push_back(obl + i * width);

	// Right col
	subnodenear.push_back(obl + 3 + i * width);

      }
      
    }
  }

  ~TravelTimeField()
  {
    delete [] nodes;
  }
  
  void construct_traveltime_field()
  {
    reset();
    
    initialize_fixed();

    while (near.get_node_count() > 0) {

      //
      // Compute trial travel times for near
      //
      update_near();

      // if (!near.validate()) {
      //  	fprintf(stderr, "error: near out of order\n");
      //  	near.dump();
      //  	throw TRAVELTIMEEXCEPTION("Near out of order\n");
      // }

      //
      // Get and fix best travel time
      //
      node_t *n = near.pop();

      n->state = node_t::STATE_FIXED;

      //
      // Add neighbours of newly fixed node to near
      //
      add_near(n);

    }
  }

  real get_traveltime(const coordinate &c)
  {
    double hnx, hny;

    coordinate::normalized_coordinate(cmin, cmax, c, hnx, hny);

    int ix = (int)((hnx - nminx)/(nmaxx - nminx) * (real)(width - 1));
    
    real alpha = (hnx - (nminx + (real)ix*dx))/dx;

    int iy = (int)((hny - nminy)/(nmaxy - nminy) * (real)(height - 1));
    real beta = (hny - (nminy + (real)iy*dy))/dy;

    real ta = (1.0 - alpha)*nodes[iy*width + ix].T + alpha*nodes[iy*width + ix + 1].T;
    real tb = (1.0 - alpha)*nodes[(iy + 1)*width + ix].T + alpha*nodes[(iy + 1)*width + ix + 1].T;

    return (1.0 - beta)*ta + beta*tb;
  }

  real get_traveltime_weighted(const coordinate &c, VelocityWeights<real> &weights)
  {
    double hnx, hny;

    coordinate::normalized_coordinate(cmin, cmax, c, hnx, hny);

    int ix = (int)((hnx - nminx)/(nmaxx - nminx) * (real)(width - 1));
    
    real alpha = (hnx - (nminx + (real)ix*dx))/dx;

    int iy = (int)((hny - nminy)/(nmaxy - nminy) * (real)(height - 1));
    real beta = (hny - (nminy + (real)iy*dy))/dy;

    weights.reset();

    
    real ta = (1.0 - alpha)*nodes[iy*width + ix].T + alpha*nodes[iy*width + ix + 1].T;
    real tb = (1.0 - alpha)*nodes[(iy + 1)*width + ix].T + alpha*nodes[(iy + 1)*width + ix + 1].T;

    nodes[iy*width + ix].back_project();
    weights.merge(nodes[iy*width + ix].vweights, (1.0 - alpha) * (1.0 - beta));
    
    nodes[iy*width + ix + 1].back_project();
    weights.merge(nodes[iy*width + ix + 1].vweights, (alpha) * (1.0 - beta));

    nodes[(iy + 1)*width + ix].back_project();
    weights.merge(nodes[(iy + 1)*width + ix].vweights, (1.0 - alpha) * (beta));

    nodes[(iy + 1)*width + ix + 1].back_project();
    weights.merge(nodes[(iy * 1)*width + ix + 1].vweights, (alpha) * (beta));

    return (1.0 - beta)*ta + beta*tb;
  }

  real get_traveltime_cubic(const coordinate &c)
  {
    double hnx, hny;

    coordinate::normalized_coordinate(cmin, cmax, c, hnx, hny);

    int ix = (int)((hnx - nminx)/(nmaxx - nminx) * (real)(width - 1));
    
    real alpha = (hnx - (nminx + (real)ix*dx))/dx;

    int iy = (int)((hny - nminy)/(nmaxy - nminy) * (real)(height - 1));
    real beta = (hny - (nminy + (real)iy*dy))/dy;

    real t0 = cerp(nodes[(iy - 1)*width + ix - 1].T,
		   nodes[(iy - 1)*width + ix].T,
		   nodes[(iy - 1)*width + ix + 1].T,
		   nodes[(iy - 1)*width + ix + 2].T,
		   alpha);

    real t1 = cerp(nodes[(iy)*width + ix - 1].T,
		   nodes[(iy)*width + ix].T,
		   nodes[(iy)*width + ix + 1].T,
		   nodes[(iy)*width + ix + 2].T,
		   alpha);
    
    real t2 = cerp(nodes[(iy + 1)*width + ix - 1].T,
		   nodes[(iy + 1)*width + ix].T,
		   nodes[(iy + 1)*width + ix + 1].T,
		   nodes[(iy + 1)*width + ix + 2].T,
		   alpha);

    real t3 = cerp(nodes[(iy + 2)*width + ix - 1].T,
		   nodes[(iy + 2)*width + ix].T,
		   nodes[(iy + 2)*width + ix + 1].T,
		   nodes[(iy + 2)*width + ix + 2].T,
		   alpha);
    
    return cerp(t0, t1, t2, t3, beta);
  }

  static real cerp(real z0, real z1, real z2, real z3, real alpha)
  {
    /* Cubic */
    real a0,a1,a2,a3;
    
    a0 = z3 - z2 - z0 + z1;
    a1 = z0 - z1 - a0;
    a2 = z2 - z0;
    a3 = z1;

    return
      a3 + alpha*(a2 + alpha*(a1 + alpha*a0));
      // a0 * alpha*alpha*alpha +
      // a1 * alpha*alpha +
      // a2 * alpha +
      // a3;
  }

  static real cerp_weights(real z0, real z1, real z2, real z3, real alpha,
			   real &w0, real &w1, real &w2, real &w3)
  {
  }
  

  bool save_traveltime_field(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }

    for (size_t j = 0; j < height; j ++) {
      for (size_t i = 0; i < width; i ++) {

	fprintf(fp, "%15.9f ", (double)nodes[j * width + i].T);

      }
      fprintf(fp, "\n");
    }

    fclose(fp);
    return true;
  }
    

  void reset()
  {
    if (refinement > 0) {
      subfield->update(velocity.field);
    }
    
    size_t size = width * height;
    for (size_t i = 0; i < size; i ++) {

      nodes[i].reset(velocity);

    }

    near.reset();
  }
  
  void initialize_fixed()
  {
    //
    // Compute the surrounding 4 nodes here
    //
    if (refinement == 0) {


      //
      // Initial fixed nodes are (col,row) .. (col + 1, row + 1) so fix them
      // now.
      //
      int bl, br, tl, tr;

      find_enclosing_nodes(hypo, bl, br, tl, tr);
      
      nodes[bl].initializeT(velocity,
			    coordinate::distance_km(hypo,
						    coordinate::from_normalized_coordinate(cmin, cmax,
											   nodes[bl].nx,
											   nodes[bl].ny)));
      nodes[br].initializeT(velocity,
			    coordinate::distance_km(hypo,
						    coordinate::from_normalized_coordinate(cmin, cmax,
											   nodes[br].nx,
											   nodes[br].ny)));
      nodes[tl].initializeT(velocity,
			    coordinate::distance_km(hypo,
						    coordinate::from_normalized_coordinate(cmin, cmax,
											   nodes[tl].nx,
											   nodes[tl].ny)));
      nodes[tr].initializeT(velocity,
			    coordinate::distance_km(hypo,
						    coordinate::from_normalized_coordinate(cmin, cmax,
											   nodes[tr].nx,
											   nodes[tr].ny)));
      
      //
      // Now add neighbours being careful not to add any of the four
      // fixed nodes and adding the same nodes each time. We use the
      // current state of the nodes to determine this.
      //

      add_near(&nodes[bl]);
      add_near(&nodes[br]);
      add_near(&nodes[tl]);
      add_near(&nodes[tr]);
      
    } else {

      subtraveltime->construct_traveltime_field();

      //
      // Fix nodes with computed child nodes
      //
      for (auto &s : subnodemappings) {
	nodes[s.first].initializeT(subtraveltime->nodes[s.second].T);
	// printf("%d %10.6f %p %d\n", s.first, nodes[s.first].T, &nodes[s.first], &nodes[s.first] - nodes);
      }

      //
      // Add near
      //
      // printf("%d near (%d)\n", (int)subnodenear.size(), (int)dirty_count);
      for (auto &i : subnodenear) {

	// printf("  %p %d\n", &nodes[i], &nodes[i] - nodes);
	add_near(&nodes[i]);
	
      }
      // printf("added (%d)\n", (int)dirty_count);
      for (size_t i = 0; i < dirty_count; i ++) {
	// printf("  %p %d\n", dirty[i], dirty[i] - nodes);
      }
    }
  }

  void add_near(node_t *fixed_node)
  {
    for (int i = 0; i < 4; i ++) {
      if (fixed_node->neighbors[i].neighbor != nullptr) {

	node_t *n = fixed_node->neighbors[i].neighbor;

	switch(n->state) {
	case node_t::STATE_FAR:
	  // Added to near
	  n->state = node_t::STATE_NEAR;
	  n->T = LARGE_T;
	  near.add_near(n);
	  mark(n);
	  break;

	case node_t::STATE_NEAR:
	  // Needs updating of travel time
	  mark(n);
	  break;

	case node_t::STATE_FIXED:
	  // Nothing
	  break;

	default:
	  // Error
	  break;
	}
      }
    }
  }
  
  void update_near()
  {
    for (size_t i = 0; i < dirty_count; i ++) {

      dirty[i]->T = dirty[i]->computeT(velocity);

      near.update_near(dirty[i]);

    }

    dirty_count = 0;
  }

  void mark(node_t *n)
  {
    //
    // Don't double add
    //
    for (size_t i = 0; i < dirty_count; i ++) {
      if (dirty[i] == n) {
	return;
      }
    }

    if (dirty_count == DIRTY_SIZE) {
      throw TRAVELTIMEEXCEPTION("Dirty array full");
    }

    dirty[dirty_count] = n;
    dirty_count ++;
  }

  void find_enclosing_nodes(const coordinate &hypo,
			    int &bl,
			    int &br,
			    int &tl,
			    int &tr)
  {
    double hnx, hny;
    
    coordinate::normalized_coordinate(cmin, cmax, hypo, hnx, hny);
    
    double cwidth = 1.0/(double)(width - 2);
    double cheight = 1.0/(double)(height - 2);
    
    int col = (int)((hnx + cwidth/2.0)/(1.0 + 2.0*cwidth) * (double)width);
    int row = (int)((hny + cheight/2.0)/(1.0 + 2.0*cheight) * (double)height);

    bl = row * width + col;
    br = row * width + col + 1;
    tl = (row + 1) * width + col;
    tr = (row + 1) * width + col + 1;
  }

  void get_traveltime_image(double *tt) const
  {
    for (size_t j = 0; j < height; j ++) {
      for (size_t i = 0; i < width; i ++) {

	tt[j * width + i] = nodes[j * width + i].T;

      }
    }
  }
	

  //
  //
  //
  coordinate cmin;
  coordinate cmax;
  coordinate hypo;

  velocity_t &velocity;
  
  size_t width;
  size_t height;

  int refinement;
  TravelTimeField *subtraveltime;
  Subfield<real> *subfield;
  velocity_t *subvelocity;
  typedef std::pair<int, int> subnodemap_t;
  std::vector<subnodemap_t> subnodemappings;
  std::vector<int> subnodenear;
  
  //
  // Near is a ordered circular buffer
  //
  CircularBuffer<real> near;
  node_t *dirty[DIRTY_SIZE];
  size_t dirty_count;

  //
  // Mesh of travel time nodes
  //
  node_t *nodes;
  real dx;
  real dy;
  real nminx;
  real nmaxx;
  real nminy;
  real nmaxy;
  
};

#endif // traveltimefield_hpp
