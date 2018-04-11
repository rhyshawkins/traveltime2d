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
#ifndef circularbuffer_hpp
#define circularbuffer_hpp

#include "traveltimenode.hpp"

#include "traveltimeexception.hpp"

template
<
  typename real
>
class CircularBuffer {
public:

  typedef TravelTimeNode<real> node_t;

  CircularBuffer(size_t maxnodes) :
    size(maxnodes),
    head(0),
    tail(0),
    full(false),
    nodes(new node_t*[maxnodes])
  {
  }

  ~CircularBuffer()
  {
    delete [] nodes;
  }

  void reset()
  {
    head = 0;
    tail = 0;
    full = false;
  }

  size_t get_size() const
  {
    return size;
  }

  size_t get_node_count() const
  {
    if (tail < head) {
      return tail + size - head;
    } else {
      if (full) {
	return size;
      } else {
	return tail - head;
      }
    }
  }
  
  node_t *pop()
  {
    if (!full && head == tail) {
      return nullptr;
    } else {
      node_t *r = nodes[head];
      head = (head + 1) % size;
      full = false;
      return r;
    }
  }

  void add_near(node_t *n)
  {
    if (n == nullptr) {
      throw TRAVELTIMEEXCEPTION("Null node");
    }
    
    if (full) {
      throw TRAVELTIMEEXCEPTION("Full circular buffer");
    }
    
    nodes[tail] = n;
    n->set_order(tail);
    tail = (tail + 1) % size;
    if (tail == head) {
      full = true;
    }
  }

  void update_near(node_t *n)
  {
    size_t i = n->get_order();

    if (nodes[i] != n) {
      throw TRAVELTIMEEXCEPTION("Order/index mismatch");
    }

    if (nodes[i]->T < 0.0) {
      throw TRAVELTIMEEXCEPTION("Negative time (%d, %f)", (int)i, nodes[i]->T);
    }

    if (!demote(i)) {
      promote(i);
    }
  }

  bool promote(size_t i)
  {
    size_t p;
    size_t c;
    
    p = (i + size - 1) % size;
    c = 0;
    while (i != head && nodes[p]->T > nodes[i]->T) {

      exchange(i, p);
      i = p;
      p = (i + size - 1) % size;
      c ++;
      
    }

    return c > 0;
  }

  bool demote(size_t i)
  {
    size_t p;
    size_t c;
    
    p = (i + 1) % size;
    c = 0;
    while (p != tail && nodes[p]->T < nodes[i]->T) {

      exchange(i, p);
      i = p;
      p = (i + 1) % size;
      c ++;
      
    }

    return c > 0;
  }

  bool validate()
  {
    size_t i = head;
    size_t next = (i + 1) % size;

    if (head == tail && !full) {
      return true;
    }

    while (next != tail) {

      if (nodes[i]->T > nodes[next]->T) {
	printf("%f %f\n", nodes[i]->T, nodes[next]->T);
	return false;
      }

      i = next;
      next = (i + 1) % size;
    } 

    return true;
  }

  void dump()
  {
    size_t i = head;
    if (head == tail && !full) {
      return;
    }

    do {
      printf("%3d %f\n", (int)i, nodes[i]->T);
      i = (i + 1) % size;
    } while (i != tail);
  }
  
private:

  void exchange(size_t i, size_t j)
  {
    node_t *t = nodes[i];
    nodes[i] = nodes[j];
    nodes[j] = t;

    nodes[i]->set_order(i);
    nodes[j]->set_order(j);
  }

  size_t size;
  size_t head;
  size_t tail;
  bool full;

  node_t **nodes;
};

#endif // circularbuffer_hpp

  
