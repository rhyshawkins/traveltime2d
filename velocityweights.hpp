#pragma once

#ifndef velocityweights_hpp
#define velocityweights_hpp

template
<
  typename real
> class VelocityWeights {
public:
  
  VelocityWeights() :
    s(16),
    n(0),
    indices(new int[16]),
    weights(new real[16]),
    dirty(true)
  {
  }
  
  ~VelocityWeights()
  {
    delete [] indices;
    delete [] weights;
  }
  
  void reset()
  {
    n = 0;
    dirty = true;
  }
  
  void add(int idx, real w)
  {
    for (int i = 0; i < n; i ++) {
      if (indices[i] == idx) {
	weights[i] += w;
	return;
      }
    }
    
    if (n == s) {
      int news = s*2;
      int *newindices = new int[news];
      real *newweights = new real[news];
      
      for (int i = 0; i < n; i ++) {
	newindices[i] = indices[i];
	newweights[i] = weights[i];
      }
      
      delete [] indices;
      delete [] weights;
      
      indices = newindices;
      weights = weights;
      
      s = news;
    }
      
    indices[n] = idx;
    weights[n] = w;
    n ++;
  }
  
  void merge(const struct VelocityWeights &vw,
	     real w)
  {
    for (int i = 0; i < vw.n; i ++) {
      add(vw.indices[i], vw.weights[i] * w);
    }
  }
  
  int s;
  int n;
  int *indices;
  real *weights;
  
  bool dirty;
};
  
#endif // velocityweights_hpp
