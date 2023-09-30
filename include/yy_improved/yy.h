//Yinyang method for accelerating exact kmeans

#ifndef YYIMP
#define YYIMP

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
#include "include/yy_improved/yy_compute_centers.h"
#include "lsh.h"

template<typename T>
struct YinyangImproved {

  typedef YyStructsImp<T> ys; //for ease of use

  //initialize the groups by running a naive kmeans a few times
  void init_groups(float* c, size_t k, size_t d, size_t t, 
    Distance& D, parlay::sequence<typename ys::center>& centers,
    parlay::sequence<typename ys::group>& groups) {

    //cluster on the groups initially using NaiveKmeans
    float* group_centers = new float[t * d];
    size_t* group_asg = new size_t[k];

    LSH<float> init;
    init(c,k,d,t,group_centers,group_asg,D);

    kmeans_bench logger = kmeans_bench(k,d,t,5,0,"LSH", "Internal Naive");
    logger.start_time();
    
    
    NaiveKmeans<float> run;
    run.cluster(c,k,d,t,
    group_centers, group_asg,D,logger, 5, 0.0001,true);
    logger.end_time();

    parlay::parallel_for(0,k,[&] (size_t i) {
      centers[i].group_id = group_asg[i];
    });

    //set group ids
    for (size_t i = 0; i < t; i++) {
      groups[i].id = i;
    }

    //TODO do in parallel?
    //sequential assigning the groups their centers
    for (size_t i = 0; i < k; i++) {
      groups[centers[i].group_id].center_ids.push_back(i);
      
    }

    delete[] group_centers; //memory cleanup
    delete[] group_asg;


  }

  //confirm that the groups are nonempty
  void assert_proper_group_size(parlay::sequence<typename ys::group>& groups, 
   size_t t, bool DEBUG_FLAG=false) {
    
    for (size_t i =0 ;i < t; i++) {
      if (groups[i].center_ids.size() == 0) {
        std::cout << 
        "Group assignment went wrong, group is wrong size"
        << std::endl;
        std::cout << groups[i].center_ids.size() << std::endl;
        abort();
      }
    }
  }

  //computes the drift, then initializes the new centers
  float update_centers_drift(size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, Distance& D, 
  parlay::sequence<typename ys::group>& groups, size_t t, float* center_calc_float, parlay::sequence<float>& deltas) {
   
    // Check convergence
    parlay::parallel_for (0,k,[&] (size_t i) { 
      centers[i].delta = sqrt_dist(
        parlay::make_slice(centers[i].coordinates).begin(), 
      center_calc_float+i*d,d,D);
      deltas[i] = centers[i].delta;
    });
    //max_diff is the largest center movement
    float max_diff = *parlay::max_element(deltas);

    //Copy over new centers
    parlay::parallel_for(0,k*d,[&] (size_t i) {
      size_t j = i/d;
      size_t coord = i%d;
      centers[j].coordinates[coord] = center_calc_float[i];
    });

    //for each group, get max drift for group
    parlay::parallel_for(0,t,[&] (size_t i) {
      auto drifts = parlay::map(groups[i].center_ids, [&] (size_t j) {
      return centers[j].delta; });
      
      groups[i].max_drift = *max_element(drifts);

    });

    return max_diff;

  }

  //update the lb's of a point given the drifts of each group
  void set_point_global_lb(typename ys::point& p, parlay::sequence<typename ys::group>& groups,
  size_t t) {
    p.global_lb = std::numeric_limits<float>::max();
    for (size_t j = 0; j < t; j++) {
      //need to cast to float for max to work
      p.lb[j] = std::max(static_cast<float>(0), p.lb[j]-groups[j].max_drift);
      //reduce the global lower bound if possible
      p.global_lb = std::min(p.global_lb,p.lb[j]);
   
    }

  }

  //yinyang needs sqrt dist for triangle inequality to work
  float sqrt_dist(float* a, float* b, size_t d, Distance& D) {
    return std::sqrt(D.distance(a,b,d));
  }

  void assert_members_n(size_t n, size_t k, const parlay::sequence<typename ys::center>& centers) {
    size_t elt_counter = 0;
    for (size_t i = 0; i < k; i++) {
      elt_counter += centers[i].old_num_members;

    }
    if (elt_counter != n) {
      std::cout << "error in num_members assignment: " << elt_counter << std::endl;
      abort();
    }
    std::cout << "Passed elt counter test " << std::endl;
  }

  //run yy
  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon,bool suppress_logging=false) {

    parlay::internal::timer tim = parlay::internal::timer();
    tim.start();
    float assignment_time = 0;
    float update_time = 0;
    float setup_time = 0;

    std::cout << "Thus begins yinyang debugging" << std::endl;

    //when we do no iterations, nothing needs to happen at all
    if (max_iter == 0) {
      return;
    }

    //we copy a point into a float buffer of fixed size 2048, so we don't support d > 2048
    if (d > 2048) {
      std::cout << "d greater than 2048, too big, printing d: " << 
      d << std::endl;
      abort();
    }
    
    //create the centers
    parlay::sequence<typename ys::center> centers = parlay::tabulate<typename ys::center>(k, [&] (size_t i) {
    return typename ys::center(i, parlay::sequence<float>(d));
    });

    //fill in the centers
    parlay::parallel_for(0,k, [&] (size_t i) {
      for (size_t j = 0; j < d; j++) {
        centers[i].coordinates[j] = *(c + i*d + j);
      }
    });
    
    //We want t to be big without overloading memory
    //because our memory cost contains O(nt) 
    size_t t;
    if (k <= 100) {t=k;}    
    else if (k <= 500) {t = k/10;}  
    else if (k<= 5000) {t = k/20;}
    else {t = k/100;}
  
    std::cout << "t is " << t << std::endl;

    
    //initialize the groups
    parlay::sequence<typename ys::group> groups(t);
    init_groups(c,k,d,t,D,centers,groups);
    assert_proper_group_size(groups,t,false); //confirm groups all nonempty

    //Init the points and point bounds:
    parlay::sequence<typename ys::point> pts = parlay::tabulate<typename ys::point>(n, [&] (size_t i) {
      return typename ys::point(asg[i],parlay::slice(v+i*d, v+i*d + d));

    });

    parlay::parallel_for(0,n,[&] (size_t i) {
  
      //first, we find the closest center to each point
      float buf[2048];
      T* it = pts[i].coordinates.begin();
      for (size_t j = 0; j < d; j++) {
          buf[j]=*it;
          it += 1; //increment the pointer
      }
      //TODO (IDEA) IMP: use an approximate square root method instead of a square root? 
      auto distances = parlay::delayed::map(centers, [&](typename ys::center& q) {
          return std::sqrt(D.distance(buf, make_slice(q.coordinates).begin(),d));
      });

      pts[i].best = min_element(distances) - distances.begin();
      pts[i].old_best = pts[i].best;

      pts[i].ub = distances[pts[i].best];

      pts[i].lb = parlay::sequence<float>(t);
      pts[i].lb[centers[pts[i].best].group_id] = std::numeric_limits<float>::max();

      for (size_t j = 0; j < k; j++) {
        if (j != pts[i].best) {
          pts[i].lb[centers[j].group_id] = 
          std::min(pts[i].lb[centers[j].group_id],distances[j]);
        }

      }
    });

    assignment_time = tim.next_time();

    //compute num_members for each center
    auto rang = parlay::delayed_tabulate(n,[] (size_t i) {return i;});

    auto new_centers_dist0 = histogram_by_key(parlay::map(rang,[&] (int i) {
      return pts[i].best;
    }));

    parlay::parallel_for(0,k,[&] (size_t i) {
      
      centers[new_centers_dist0[i].first].new_num_members=new_centers_dist0[i].second;
      centers[new_centers_dist0[i].first].old_num_members=new_centers_dist0[i].second;
      //on first iter, all centers have changed
      centers[new_centers_dist0[i].first].has_changed = true;

    });
    //debugging (confirm all points belong to a center)
    assert_members_n(n,k,centers);
    
    //iters start at 1 as we have already done a closest point check
    size_t iters = 1; 
    float max_diff = 0.0;
    //keep track of the number of distance calculations
    parlay::sequence<size_t> distance_calculations(n,k); 
    //keep track of the number of points reassigned in an iteration
    parlay::sequence<uint8_t> center_reassignments(k,1);
    //wasteful to keep track of deltas here AND in centers? TODO
    //currently, have this array because easier to manipulate/do operations on (instead of doing a parlay::map from centers)
    //which is more efficient?
    parlay::sequence<float> deltas(k);
    //Step 3: Repeat until convergence

    //for center calculation
    float* center_calc_float = new float[k*d];
    //TODO remove center_calc?
    //TODO stop static casting to double?
    double* center_calc = new double[k*d];

    //in the first iteration, we use a standard compute centers
    //whereas in future iterations, we theoretically would
    //use a comparative compute centers
    bool first_time = true;
    
    //Compute Centers Object
    YyComputeCentersImp<T> comp_cen;

    setup_time = tim.next_time();

    //our iteration loop, will stop when we've done max_iter iters, or if we converge (within epsilon)
    while (true) {

      //the first iteration, we don't have previous centers to compare to
      //so we must do a full-adding version of compute centers
      if (first_time) {

         comp_cen.compute_centers_filter(v,pts,n,d,k,centers,center_calc,center_calc_float);
         
        first_time=true; //TODO set to false to use the comparative compute centers method
      }
      else {
        comp_cen.compute_centers_comparative5(pts,n,d,k,centers,center_calc,center_calc_float);
      }

      max_diff = update_centers_drift(d,k,centers,D,groups,t,center_calc_float,deltas);
     
      //TODO is this correct? (for getting the msse's?)
      //because the ub may not be tight
      parlay::sequence<double> best_dists = parlay::tabulate(n,[&] (size_t i) {
        return static_cast<double>(pts[i].ub * pts[i].ub);
      });
      double squared_error = parlay::reduce(best_dists)/n; //divide by n for msse

      update_time = tim.next_time();

      //end of iteration stat updating
      if (!suppress_logging) {
         logger.add_iteration(assignment_time,update_time,squared_error,parlay::reduce(distance_calculations),
      parlay::reduce(center_reassignments),deltas,setup_time);


      }
      assignment_time=0;
      update_time=0;
      setup_time=0;
      //convergence check
      if (iters >= max_iter || max_diff <= epsilon) break;

      iters += 1;
      distance_calculations = parlay::sequence<size_t>(n,0); //reset distance calc counter
      center_reassignments = parlay::sequence<uint8_t>(n,0);
      if (!suppress_logging) {
        std::cout << "iter: " << iters << std::endl;
      }

      //set centers changed to false, update old_num_members
      parlay::parallel_for(0,k,[&] (size_t i) {
        centers[i].has_changed = false;
        centers[i].old_num_members=centers[i].new_num_members;
      });

      //3.2: Group filtering      
      parlay::parallel_for(0,n,[&](size_t i) {

        //update bounds and old_best
        pts[i].ub += centers[pts[i].best].delta; 
        pts[i].old_best = pts[i].best; 
        set_point_global_lb(pts[i],groups,t);

        //nothing happens if our closest center can't change
        //let's tighten our ub every iteration (by not doing an initial check here)
        // if (pts[i].global_lb >= pts[i].ub) {
        //   return;
        // }
      

        //copy to float buffer
        float buf[2048];
        T* it = pts[i].coordinates.begin();
        for (size_t coord = 0; coord < d; coord++) {
          buf[coord] = *it;
          it += 1;
        }

        //tighten the upper bound
        pts[i].ub = sqrt_dist(buf,
        parlay::make_slice(centers[pts[i].best].coordinates).begin(),
        d,D);
        distance_calculations[i] += 1;

        //again, nothing happens if our closest center can't change
        if (pts[i].global_lb >= pts[i].ub) {
          return;
        }

        //for each group
        for (size_t j = 0; j < t; j++) {
          //if group j is too far away we don't look at it
          if (pts[i].ub <= pts[i].lb[j]) {
            continue;
          }
                   
          //reset the lower bound, make it smallest distance we calculate that's not the closest distance away
          pts[i].lb[j] = std::numeric_limits<float>::max(); 
            
          //for each group member (center)
          for (size_t l = 0; l < groups[j].center_ids.size(); l++) {

              //don't do a distance comparison with the previous best
              if (pts[i].old_best == groups[j].center_ids[l]) {
                continue;
              }

              //find distance to center l in group j
              float new_d = sqrt_dist(
              buf,
              parlay::make_slice(
                centers[groups[j].center_ids[l]].coordinates).begin(),
              d,D);

              //increment distance calc counter for this pt
              distance_calculations[i]++; 

              //note that the ub is tight rn,
              //that ub IS the distance to the previously 
              //closest center
              //So if our new dist is less than ub, we have
              //a new closest point
              if (pts[i].ub > new_d) {
                //the group with the previous center gets a slightly
                //lower bound, because the distance to the old center can 
                //become a lower bound
                //minus needed because of the center change from this iter
                //TODO is this a correct adjustment of the lower bound?? consider again
                pts[i].lb[centers[pts[i].best].group_id]=
                std::max(pts[i].ub-centers[pts[i].best].delta, std::min(pts[i].ub,
                pts[i].lb[centers[pts[i].best].group_id] - centers[pts[i].best].delta));
                
                pts[i].best=groups[j].center_ids[l];
                //new ub is tight
                pts[i].ub = new_d;

                //log center reassign
                center_reassignments[i] = 1;
                //mark centers have changed. yes this is a race, but because we are setting false to true this is fine
                centers[pts[i].best].has_changed = true;
                centers[pts[i].old_best].has_changed=true;
                  
              }
              else {
                //if this center is not the closest, use it to improve the lower bound
                pts[i].lb[j] = std::min(new_d,pts[i].lb[j]);
              
              }
          }  
        }
      }); //gran 1? I think helps.

      assignment_time = tim.next_time();

      //record num_new_members for each center
      auto new_centers_dist = histogram_by_key(parlay::map(pts,[&] (typename ys::point& p) {
        return p.best;
      }));
      for (size_t i = 0; i < k; i++) {
        centers[new_centers_dist[i].first].new_num_members = new_centers_dist[i].second; 
      }
      
      setup_time = tim.next_time();

    }
   
    //copy back over coordinates
    //put our data back 
    parlay::parallel_for(0,k,[&] (size_t i) {
      for (size_t j = 0; j < d; j++) {
        c[i*d + j] = centers[i].coordinates[j];
      }
    });
    parlay::parallel_for(0,n,[&] (size_t i) {
        asg[i] = pts[i].best;
    });


    delete[] center_calc;
    delete[] center_calc_float;

  }

};

#endif //YYIMP