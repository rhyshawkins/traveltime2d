#pragma once

#ifndef velocityweights_hpp
#define velocityweights_hpp

#define VW_USE_ORDERED

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
  
#ifdef VW_USE_ORDERED
  void add(int idx, real w)
  {
    if (n == 0 || indices[n - 1] < idx) {
      //
      // Append
      //
      resize(n + 1);
      indices[n] = idx;
      weights[n] = w;

      n ++;
    } else {

      for (int i = 0; i < n; i ++) {

	if (indices[i] == idx) {
	  //
	  // Add
	  //
	  weights[i] += w;
	  return;
	}
	
	if (indices[i] > idx) {
	  //
	  // Insert
	  //
	  resize(n + 1);
	  for (int j = n; j > i; j --) {
	    indices[j] = indices[j - 1];
	    weights[j] = weights[j - 1];
	  }
	  indices[i] = idx;
	  weights[i] = w;
	  n ++;
	  
	  return;
	}
      }

      printf("Failed to add index\n");
    }
  }
#else
  void add(int idx, real w)
  {
    for (int i = 0; i < n; i ++) {
      if (indices[i] == idx) {
	weights[i] += w;
	return;
      }
    }

    resize(n + 1);
      
    indices[n] = idx;
    weights[n] = w;
    n ++;
  }
#endif // VW_USE_ORDERED

#ifdef VW_USE_ORDERED
  void merge(const struct VelocityWeights &vw,
	     real w)
  {
    int i = 0;
    int j = 0;

    while (j < vw.n) {

      if (i >= n) {
	//
	// Append
	//
	resize(n + 1);
	indices[n] = vw.indices[j];
	weights[n] = w * vw.weights[j];
	n ++;
	i ++;
	j ++;
      } else if (indices[i] == vw.indices[j]) {
	weights[i] += w * vw.weights[j];
	i ++;
	j ++;
      } else if (indices[i] < vw.indices[j]) {
	i ++;
      } else {
	//
	// Insert
	//
	resize(n + 1);
	for (int k = n; k > i; k --) {
	  indices[k] = indices[k - 1];
	  weights[k] = weights[k - 1];
	}
	n ++;
	indices[i] = vw.indices[j];
	weights[i] = w * vw.weights[j];
	i ++;
	j ++;
      }
    }
  }
  
#else
  
  void merge(const VelocityWeights &vw,
	     real w)
   {
     for (int i = 0; i < vw.n; i ++) {
       add(vw.indices[i], vw.weights[i] * w);
     }
   }
#endif // VW_USE_ORDERED

  void assign(const VelocityWeights &vw)
  {
    resize(vw.n);

    n = vw.n;
    for (int i = 0; i < n; i ++) {
      indices[i] = vw.indices[i];
      weights[i] = vw.weights[i];
    }
  }
  
  // {
  //   int i = 0;
  //   int j = 0;

  //   while (j < vw.n) {

  //     if (i >= n) {
  // 	//
  // 	// Append
  // 	//

  // 	j ++;
  //     } else if (indices[i] == vw.indices[j]) {

  // 	weights[i] += vw.weights[j];
  // 	i ++;
  // 	j ++;

  //     } else if (indices[i] < vw.indices[j]) {

  // 	i ++;

  //     } else {
  // 	//
  // 	// Insert
  // 	//
	
  //     }

  //   }
  // }

  // int insertion_index(int idx,
  // 		      int start,
  // 		      int end)
  // {
  //   if (start > end) {
  //     return start;
  //   }

  //   if (idx <= indices[start]) {
  //     return start;
  //   }

  // }

  // void insert_at(int i,
  // 		 int idx,
  // 		 double weight)
  // {
  //   if (i < n && indices[i] == idx) {
  //     weights[i] += weight;
  //     return;
  //   }

  //   if (n == s) {
  //     //
  //     // resize
  //     //
  //     resize(n + 1);
  //   }
      
  //   for (int j = n; j > i; j --) {
  //     indices[j] = indices[j - 1];
  //     weights[j] = weights[j - 1];
  //   }
  // }

  //   void resize(int requied)
  //   {
  //     while (s < required) {
  //     }
  //   }

  void resize(int required)
  {
    int news = s;

    while (news < required) {
      news *= 2;
    }

    if (news != s) {
      int *newindices = new int[news];
      real *newweights = new real[news];
      
      for (int i = 0; i < n; i ++) {
	newindices[i] = indices[i];
	newweights[i] = weights[i];
      }
      
      delete [] indices;
      delete [] weights;
      
      indices = newindices;
      weights = newweights;
      
      s = news;
    }
  }
  
  int s;
  int n;
  int *indices;
  real *weights;
  
  bool dirty;
};
  
#endif // velocityweights_hpp
