//compute centers for yy getting messy, so
//moving to its own file
#ifndef YYCENTERS
#define YYCENTERS


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
#include "include/yy_structs.h"


template <typename T>
struct YyComputeCenters {


  typedef YyStructs<T> ys; //for ease of use

  //compute centers calculates the new centers
  //puts the new values into the float* center_calc
  //TODO make a different run (no casting) if T is a float
  //typename needed because point a dependent type
  static void compute_centers_filter(
  parlay::sequence<typename YyStructs<T>::point>& pts, size_t n, size_t d, size_t k, 
  const parlay::sequence<typename YyStructs<T>::center>& centers, double* center_calc,
  float* center_calc_float,bool suppress_logging=false) {

    //update centers timer
    //commenting out timer prints because too much printing
    parlay::internal::timer t2;
    t2.start();
    if (!suppress_logging) {
      std::cout << "using the filter version " << std::endl;
    }

    parlay::sequence<parlay::sequence<size_t>> indices(k);

    //t2.next("Made new centers");

    parlay::parallel_for(0,k,[&] (size_t i) {
      auto temp = parlay::filter(pts,[&] (typename YyStructs<T>::point& p) {
        return p.best==i;
      });
      indices[i] = parlay::map(temp, [&] (typename YyStructs<T>::point& p) {
        return p.id;
      });
    
    });

    //t2.next("Got indices");

    parlay::parallel_for (0, k*d, [&] (size_t icoord){
    size_t i = icoord / d;
    size_t coord = icoord % d;

    //if there are no values in a certain center, just keep the center 
    //where it is
    if (indices[i].size() > 0) { 
    //note the static cast to double here to avoid overflow
    center_calc[icoord] = static_cast<float>(reduce(parlay::map(indices[i],[&] 
    (size_t ind) {return static_cast<double>(
    pts[ind].coordinates[coord]);})) / indices[i].size()); 

    }
    else { 
    center_calc[icoord] = centers[i].coordinates[coord];
    }

    }); 

  

    //t2.next("Done with adding");

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    //t2.next("Done with copying");
  }

  // Compute centers by comparing to the previous results
  static void compute_centers_comparative(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, 
  double* center_calc,float* center_calc_float) {

    parlay::internal::timer t2;
    t2.start();
    // Compute new centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      center_calc[i] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

    });
     // t2.next("Multiplied");


    parlay::sequence<typename ys::point> changed_points = parlay::filter(pts,[&] (typename ys::point& p) {
      return p.best != p.old_best;
    });

   // t2.next("Filtered");

    parlay::sequence<parlay::sequence<typename ys::point>> add_these_all(k);
    parlay::sequence<parlay::sequence<typename ys::point>> sub_these_all(k);
    parlay::parallel_for (0,k,[&] (size_t i) {
      add_these_all[i] = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.best==i;
      });
      sub_these_all[i] = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.old_best == i;
      });

    });

    //t2.next("Extra filtering");

    parlay::parallel_for(0,k*d,[&] (size_t jcoord) {
      size_t j = jcoord / d;
      size_t coord = jcoord % d;

        parlay::sequence<typename ys::point> add_these = add_these_all[j];
        parlay::sequence<typename ys::point> sub_these = sub_these_all[j];
        double diff_to_add = 0;
        for (size_t m = 0; m < add_these.size(); m++) {
          diff_to_add += add_these[m].coordinates[coord];
        }
        for (size_t m = 0; m < sub_these.size(); m++) {
          diff_to_add -= sub_these[m].coordinates[coord];
        }
        center_calc[j*d+coord] += diff_to_add;

    }); 

      //  t2.next("Added/subbed changes");


      parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      if (centers[j].new_num_members != 0) {
        center_calc[i] /= centers[j].new_num_members;
      }
      else {
        center_calc[i] = centers[j].coordinates[dim];
      }
    });

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    //t2.next("Copied back");

  }

  // Compute centers by comparing to the previous results
  //this version will parallelise differently
  //takes .3-.4s/iter on k=1mil, n=10, iter=5
  static void compute_centers_comparative2(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, 
  double* center_calc,float* center_calc_float) {

    parlay::internal::timer t2;
    t2.start();
    // Compute new centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      center_calc[i] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

    });
    //  t2.next("Multiplied2");


    parlay::sequence<typename ys::point> changed_points = parlay::filter(pts,[&] (typename ys::point& p) {
      return p.best != p.old_best;
    });

    //t2.next("Filtered2");

    parlay::sequence<parlay::sequence<typename ys::point>> add_these_all(k);
    parlay::sequence<parlay::sequence<typename ys::point>> sub_these_all(k);
    parlay::parallel_for (0,k,[&] (size_t i) {
      add_these_all[i] = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.best==i;
      });
      sub_these_all[i] = parlay::filter(changed_points, [&] (typename ys::point& p) {
        return p.old_best == i;
      });

    });

   // t2.next("Extra filtering2");

    parlay::parallel_for(0,k*d,[&] (size_t jcoord) {
      size_t j = jcoord / d;
      size_t coord = jcoord % d;
 
    
        double diff_to_add = 0;
        for (size_t m = 0; m < add_these_all[j].size(); m++) {
          diff_to_add += add_these_all[j][m].coordinates[coord];
        }
        for (size_t m = 0; m < sub_these_all[j].size(); m++) {
          diff_to_add -= sub_these_all[j][m].coordinates[coord];
        }
        center_calc[j*d+coord] += diff_to_add;

    });

   // t2.next("Added/subbed changes");

    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      if (centers[j].new_num_members != 0) {
        center_calc[i] /= centers[j].new_num_members;
      }
    });

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    t2.next("Copied back");

  }

  // Compute centers by comparing to the previous results
  // .5 s
  void compute_centers_comparative3(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, 
  double* center_calc,float* center_calc_float) {

    parlay::internal::timer t2;
    t2.start();
    // Compute new centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      center_calc[i] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

    });
      t2.next("Multiplied3");


    parlay::sequence<typename ys::point> changed_points = parlay::filter(pts,[&] (typename ys::point& p) {
      return p.best != p.old_best;
    });

    t2.next("Filtered");

    parlay::parallel_for(0,k,[&] (size_t j) {
    
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
    },5); 

    t2.next("Added/subbed changes");

    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      if (centers[j].new_num_members != 0) {
        center_calc[i] /= centers[j].new_num_members;
      }
    });

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    t2.next("Copied back");


    }


  // Compute centers by comparing to the previous results
  // Version 4 almost fully sequential 
  // This honestly the fastest so far (.25s)? (depends on input)
  void compute_centers_comparative4(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, 
  double* center_calc,float* center_calc_float) {

    parlay::internal::timer t2;
    t2.start();
    // Compute new centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      center_calc[i] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

    });
      t2.next("Multiplied");


    parlay::sequence<typename ys::point> changed_points = parlay::filter(pts,[&] (typename ys::point& p) {
      return p.best != p.old_best;
    });

    t2.next("Filtered");


    //Remark: changed_points usually very small so this is okay
    for (size_t i = 0; i < changed_points.size(); i++) {
      //TODO should this inner loop be ||? (on d)
      typename ys::point p = pts[i];
      for (size_t coord =0 ; coord < d; coord++) {
      //TODO CAST TO DOUBLE
        center_calc[p.best * d + coord] += static_cast<double>(p.coordinates[coord]);
        center_calc[p.old_best *d + coord] -= static_cast<double>(p.coordinates[coord]);
      }

    }

    t2.next("Added/subbed changes");


      parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      if (centers[j].new_num_members != 0) {
        center_calc[i] /= centers[j].new_num_members;
      }
    });

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    t2.next("Copied back");

  }
    //compute centers 5

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
      // else {
      //   for (size_t coord = 0; coord < d; coord++) {
      //     center_calc[j*d+coord] = static_cast<double>(centers[j].coordinates[coord]);
      //   }
      // }
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


};


#endif // YYCENTERS