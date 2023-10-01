//compute centers for yy getting messy, so
//moving to its own file
#ifndef YYCENTERS_IMP
#define YYCENTERS_IMP


#include <cmath>
#include <tuple>
#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/delayed.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"
#include "utils/NSGDist.h"
#include "initialization.h"
#include "naive.h"
#include "include/utils/kmeans_bench.h"
#include "include/yy_improved/yy_structs.h"


template <typename T>
struct YyComputeCentersImp {


  typedef YyStructsImp<T> ys; //for ease of use

  

  //compute centers by subtracting/adding points that left/joined a center, instead of recalculating from scratch
  //TODO review this code
   void compute_centers_comparative5(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, 
  double* center_calc,float* center_calc_float) {

    parlay::internal::timer t2;
    t2.start();
    // Compute new centers
    parlay::parallel_for(0,k,[&](size_t j) {
      if (centers[j].has_changed) {
        for (size_t dim = 0; dim < d; dim++) {
          center_calc[j*d+dim] = static_cast<double>(centers[j].coordinates[dim]) * static_cast<double>(centers[j].old_num_members);

        }

      }
      
      
    });
    //  t2.next("Multiplied3");


    parlay::sequence<typename ys::point> changed_points = parlay::filter(pts,[&] (typename ys::point& p) {
      return p.best != p.old_best;
    });

   // t2.next("Filtered");
    //TODO this seems a quite inefficient way to 
    //parallelize the yy adding structure, think about again
    parlay::parallel_for(0,k,[&] (size_t j) {
      if (centers[j].has_changed) {
        parlay::sequence<typename ys::point> add_these = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.best==j;
      });
      parlay::sequence<typename ys::point> sub_these = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.old_best == j;
      });
     
      for (size_t coord = 0; coord < d; coord++) {
        center_calc[j*d+coord] = center_calc[j*d+coord] + 
        parlay::reduce(parlay::map(add_these,[&] (typename ys::point& p) {return static_cast<double>(p.coordinates[coord]);} )) - 
        parlay::reduce(parlay::map(sub_these, [&] (typename ys::point& p) {return static_cast<double>(p.coordinates[coord]);}));

      }

      }
    
      
    }); 

   // t2.next("Added/subbed changes");

    parlay::parallel_for(0,k,[&](size_t j) {
      
      if (centers[j].new_num_members != 0 && centers[j].has_changed) {
        for (size_t coord = 0; coord < d; coord++) {
          center_calc[j*d+coord] /= static_cast<double>(centers[j].new_num_members);
        }
      }
      
    });

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      size_t j = i/d;
      if (centers[j].new_num_members==0 || !centers[j].has_changed) {
        center_calc_float[i] = centers[j].coordinates[i%d];
      }
      else {
        center_calc_float[i] = center_calc[i];
      }
    });

    //t2.next("Copied back");


    }

  //compute centers calculates the new centers the standard way
  //puts the new values into the float* center_calc_float
  //TODO make a different run (no casting) if T is a float
  void compute_centers_filter(T* v,
  const parlay::sequence<typename ys::point>& pts, size_t n, size_t d, size_t k, 
  const parlay::sequence<typename ys::center>& centers, double* center_calc,
  float* center_calc_float,bool suppress_logging=false) {
    //copy center coords into center_calc_float
    //TODO fixme this is sus?
    parlay::parallel_for(0,k,[&] (size_t i) {
      for (size_t j = 0; j < d; j++) {
        center_calc_float[i*d+j] = 0; //not c** TODO FIXME

      }
    });

    auto rangn = parlay::delayed_tabulate(n,[&] (size_t i) { return i; });
    //TODO pass in the pts_grouped_by_center <- avoid doing a group_by twice**
    parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(pts[i].best,i);
    }));

    get_groups1(v,n,d,k,center_calc_float,pts,pts_grouped_by_center);

    add_points2(v,n,d,k,center_calc_float,pts,pts_grouped_by_center);
    //if a center loses its members, do not zero out the center,
    //just keep it where it is
    parlay::parallel_for(0,k,[&] (size_t i) {

      parlay::parallel_for(0,d,[&] (size_t coord) {
        if (pts_grouped_by_center[i].second.size() == 0) {
          
          center_calc_float[pts_grouped_by_center[i].first*d+coord] = centers[pts_grouped_by_center[i].first].coordinates[coord];

        }
    
      });
    
    });

  }
  //first part of update centers is to get the partitioning of the points (we do with a group by)
  void get_groups1(T* v, size_t n, size_t d, size_t k, float* c, const parlay::sequence<typename ys::point>& pts, parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& pts_grouped_by_center) {


    auto rangn = parlay::delayed_tabulate(n,[&] (size_t i) { return i; });
    pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
      return std::pair(pts[i].best,i);
    }));
  }

  //second part of update centers is to do the actual adding
  void add_points2(T* v, size_t n, size_t d, size_t k, float* c, const parlay::sequence<typename ys::point>& pts, parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& pts_grouped_by_center) {

    parlay::parallel_for(0,k,[&] (size_t i) {
        size_t picked_center_d = pts_grouped_by_center[i].first*d;
        for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
          size_t point_coord = pts_grouped_by_center[i].second[j]*d;
          for (size_t coord = 0; coord < d; coord++) {
            c[picked_center_d + coord] += static_cast<float>(v[point_coord + coord]);
          }
        }
        

    },1);

      parlay::parallel_for(0,k,[&] (size_t i) {

      parlay::parallel_for(0,d,[&] (size_t coord) {
        if (pts_grouped_by_center[i].second.size() > 0) {
          
          c[pts_grouped_by_center[i].first*d+coord] /= pts_grouped_by_center[i].second.size();

        }
      

      });
    
    });

  }



};


#endif // YYCENTERS_IMP