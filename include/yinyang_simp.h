//Yinyang method for accelerating exact kmeans

#ifndef YYSIMP
#define YYSIMP

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
#include "include/yy_compute_centers.h"

template<typename T>
struct YinyangSimp {

  typedef YyStructs<T> ys; //for ease of use


  //for debugging, get distance between point and center
  float pc_dist(typename ys::point& p, typename ys::center& c, Distance& D, size_t d) {
    float buf2[2048];
    T* it2 = p.coordinates.begin();
    for (size_t coord = 0; coord < d; coord++) {
      buf2[coord] = *it2;
      it2 += 1;
    }
    return sqrt_dist(buf2,
    parlay::make_slice(c.coordinates).begin(),d,D);
    //std::cout << "dist to center 16: " << dist41 << std::endl;

  }

  //copy the input data into our centers sequence
  void fill_in_centers(float* c, parlay::sequence<typename ys::center>& centers, 
  size_t k, size_t d) {
    parlay::parallel_for(0,k, [&] (size_t i) {
    for (size_t j = 0; j < d; j++) {
    centers[i].coordinates[j] = *(c + i*d + j);
    }

    });
  }

  

  //initialize the groups by running a naive kmeans a few times
  void init_groups(float* c, size_t k, size_t d, size_t t, 
    Distance& D, parlay::sequence<typename ys::center>& centers,
    parlay::sequence<typename ys::group>& groups) {

    //std::cout << "Debugging init groups " << std::endl;

    //cluster on the groups initially using NaiveKmeans
    float* group_centers = new float[t * d];
    size_t* group_asg = new size_t[k];


    LazyStart<float> init;
    init(c,k,d,t,group_centers,group_asg,D);

    kmeans_bench logger = kmeans_bench(k,d,t,5,0,"Lazy", "Internal Naive");
    logger.start_time();
    
    
    NaiveKmeans<float> run;
    run.cluster(c,k,d,t,
    group_centers, group_asg,D,logger, 5, 0.0001,true);
    logger.end_time();
    //std::cout << "survived group init" << std::endl;

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
    if (DEBUG_FLAG) {
       std::cout << "printing out init groups" << std::endl;
      for (size_t i = 0; i < t; i++) {
        ys::print_group(groups[i]);
      }
    }
    
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
  float update_centers_drift(parlay::sequence<typename ys::point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<typename ys::center>& centers, Distance& D, 
  parlay::sequence<typename ys::group>& groups, size_t t, float* center_calc_float) {
   
    // Check convergence
    parlay::sequence<float> total_diffs_arr(k,0);
    parlay::parallel_for (0,k,[&] (size_t i) { 
      centers[i].delta = sqrt_dist(
        parlay::make_slice(centers[i].coordinates).begin(), 
      center_calc_float+i*d,d,D);
      total_diffs_arr[i] = centers[i].delta;
    });

    float total_diff = parlay::reduce(total_diffs_arr);

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

    return total_diff;

  }

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

  void debug_group_asg_centers(parlay::sequence<typename ys::center>& centers,float* group_centers, float* group_asg,
  size_t n, size_t d, size_t k, size_t t) {
    std::cout << "Couting same info as above (group asg) but outside of the function" << std::endl;
    std::cout << "Group asg: " << std::endl;
    for (size_t i = 0; i < k; i++) {
      std::cout << group_asg[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Couting group centers: " << std::endl;
    for (size_t i = 0; i < t; i++) {
      for (size_t j = 0; j < d; j++) {
        std::cout << group_centers[i*d+j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "Couting group ids" << std::endl;
    for (size_t i = 0; i < k; i++) {
      std::cout << centers[i].group_id << " ";
    }
    std::cout << std::endl;

    std::cout << "made it here1" << std::endl;
  }

  //suppress logging used in quantized method
  //seeing the yy runs of the individual blocks 
  //can be confusing, especially as they are run in parallel
  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon,bool suppress_logging=false) {

    parlay::internal::timer tim = parlay::internal::timer();
    tim.start();

    std::cout << "Thus begins yinyang debugging" << std::endl;

    //when we do no iterations, nothing needs to happen at all
    if (max_iter == 0) {
      return;
    }

    if (d > 2048) {
      std::cout << "d greater than 2048, too big, printing d: " << 
      d << std::endl;
      abort();
    }
    
    //format the data
    parlay::sequence<typename ys::point> pts = parlay::tabulate<typename ys::point>(n, [&] (size_t i) {
    return typename ys::point(asg[i],parlay::slice(v+i*d, v+i*d + d));

    });

    parlay::parallel_for(0,n,[&] (size_t i) {
      pts[i].id=  i;
    });

    //create the centers
    parlay::sequence<typename ys::center> centers = parlay::tabulate<typename ys::center>(k, [&] (size_t i) {
    return typename ys::center(i, parlay::sequence<float>(d));
    });

    //fill in the centers
    fill_in_centers(c,centers,k,d);
    
    //We want t to be big without overloading memory
    //because our memory cost contains O(nt) 
    size_t t;
    if (k <= 100) {
      t=k;
    }
    else if (k <= 1000) {
      t = k/10;

    }
    else if (k <= 10000) {
      t = k/20;
    }
    else {
      t = k/100;
    }
  
    parlay::sequence<typename ys::group> groups(t);
    //initialize the groups
    init_groups(c,k,d,t,D,centers,groups);
    if (!suppress_logging) {
       logger.add_iteration(0,0,45,0,0,
    parlay::sequence<float>(k,0),tim.next_time());
    }
   

    //confirm group assignment happened properly
    //(checking a necessary not sufficient condition)
    assert_proper_group_size(groups,t);

    //Init the point bounds:

    //first, we find the closest center to each point
    parlay::parallel_for(0,n,[&] (size_t i) {
       float buf[2048];
        T* it = pts[i].coordinates.begin();
        for (size_t j = 0; j < d; j++) {
            buf[j]=*it;
            it += 1; //increment the pointer
        }
        
        auto distances = parlay::delayed::map(centers, [&](typename ys::center& q) {
            return std::sqrt(D.distance(buf, make_slice(q.coordinates).begin(),d));
        });

        pts[i].best = min_element(distances) - distances.begin();
        pts[i].old_best = pts[i].best;

        pts[i].ub = distances[pts[i].best];
        pts[i].lb = parlay::sequence<float>(t,std::numeric_limits<float>::max());

        for (size_t j = 0; j < k; j++) {
          if (j != pts[i].best) {
            pts[i].lb[centers[j].group_id] = 
            std::min(pts[i].lb[centers[j].group_id],distances[j]);
          }

        }
    });
    if (!suppress_logging) {
       logger.add_iteration(0,0,48,0,0,
    parlay::sequence<float>(k,0),tim.next_time());

    }
   

    //set num_members for the centers
    parlay::parallel_for(0,k,[&] (size_t i) {
      auto belonging_pts = parlay::filter(pts,[&] (typename ys::point& p) {
        return p.best==i;
      });
      centers[i].new_num_members=belonging_pts.size();
      centers[i].old_num_members=belonging_pts.size();
    });
    if (!suppress_logging) {
      logger.add_iteration(0,0,49,0,0,parlay::sequence<float>(k,0),
    tim.next_time());

    }
    

    //debugging
    {
    size_t elt_counter = 0;
    for (size_t i = 0; i < k; i++) {
      elt_counter += centers[i].old_num_members;

    }
    if (elt_counter != n) {
      std::cout << "error in num_members assignment: " << elt_counter << std::endl;
      abort();
    }}
    std::cout << "Passed elt counter test " << std::endl;

    //iters start at 1 as we have already done a closest point check
    size_t iters = 1; 
    float total_diff = 0.0;
    //keep track of the number of distance calculations
    parlay::sequence<size_t> distance_calculations = parlay::sequence<size_t>(k);
    //keep track of the number of points reassigned in an iteration
    parlay::sequence<uint8_t> center_reassignments = parlay::sequence<uint8_t>(n);

    float assignment_time = tim.next_time();
    float update_time = 0;
    float setup_time = 0;
    //Step 3: Repeat until convergence

    //for center calculation
    double* center_calc = new double[k*d];
    float* center_calc_float = new float[k*d];

    bool first_time = true;
    
    //Compute Centers Object
    YyComputeCenters<T> comp_cen;

    //on the first iteration, all centers have changed
    parlay::parallel_for(0,k,[&] (size_t i) {
      centers[i].has_changed = true;
    });

    
    while (true) {
      //the first iteration, we don't have previous centers to compare to
      //so we must do a full-adding version of compute centers
      if (first_time) {
        comp_cen.compute_centers_filter(pts,n,d,k,centers,center_calc,center_calc_float,suppress_logging);
        first_time=false;
      }
      else {
        comp_cen.compute_centers_comparative5(pts,n,d,k,centers,center_calc,center_calc_float);
      }
      total_diff = update_centers_drift(pts,n,d,k,centers,D,groups,t,center_calc_float);
      
      parlay::sequence<float> deltas = parlay::tabulate(k,[&] (size_t i) {
        return centers[i].delta;
      });
      parlay::sequence<double> best_dists = parlay::tabulate(n,[&] (size_t i) {
        return static_cast<double>(pts[i].ub * pts[i].ub);
      });
      double squared_error = parlay::reduce(best_dists);

      update_time = tim.next_time();

      //end of iteration stat updating
      if (!suppress_logging) {
         logger.add_iteration(assignment_time,update_time,squared_error,parlay::reduce(distance_calculations),
      parlay::reduce(center_reassignments),deltas,setup_time);


      }
      //convergence check
      if (iters >= max_iter || total_diff <= epsilon) break;

      iters += 1;
      distance_calculations = parlay::sequence<size_t>(n,0); //reset distance calc counter
      center_reassignments = parlay::sequence<uint8_t>(n,0);
      if (!suppress_logging) {
        std::cout << "iter: " << iters << std::endl;
      }

      //on the first iteration, all centers have changed
      parlay::parallel_for(0,k,[&] (size_t i) {
        centers[i].has_changed = false;
        centers[i].old_num_members=centers[i].new_num_members;
      });

      //3.2: Group filtering

      parlay::parallel_for(0,n,[&](size_t i) {
        //update bounds
        pts[i].ub += centers[pts[i].best].delta; 
        //set old best here, as the old best is independent of 
        //whether we change
        pts[i].old_best = pts[i].best; 
        
        set_point_global_lb(pts[i],groups,t);
      
        
        //if our center could change
        if (pts[i].global_lb < pts[i].ub) {

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

          //check the condition again
          if (pts[i].global_lb < pts[i].ub) {
           
            //TODO add in local filtering 
           
            //for each group
            for (size_t j = 0; j < t; j++) {
              //if the lower bound of group j is smaller than the 
              //upper bound of pt i
              if (pts[i].ub > pts[i].lb[j]) {
                //look at each group member (ie each center)

                //reset the lower bound, make it smallest distance we calculate that's not the closest distance away
                pts[i].lb[j] = std::numeric_limits<float>::max(); 
                
                //for each group member (center)
                for (size_t k = 0; k < groups[j].center_ids.size(); k++) {

                    //since we skip the local filter, we do a distance calculation
                    float new_d = sqrt_dist(
                    buf,
                    parlay::make_slice(
                      centers[groups[j].center_ids[k]].coordinates).begin(),
                    d,D);

                    //increment distance calc counter
                    distance_calculations[i]++; 

                    //note that the ub is tight rn,
                    //that ub IS the distance to the previously 
                    //closest center
                    //So if our new dist is less than ub, we have
                    //a new closest point
                    if (pts[i].ub > new_d) {
                      
                      //mark point has changed, adjust center assignment counts
                      //shoot this isn't thread safe
                      //centers[groups[j].center_ids[k]].new_num_members += 1;
                      //centers[pts[i].best].new_num_members -= 1;

                    

                      //the group with the previous center gets a slightly
                      //lower bound, because the distance to the old center can 
                      //become a lower bound
                      //minus needed because of the center change from this iter
                      pts[i].lb[centers[pts[i].best].group_id]=
                      std::max(pts[i].ub-centers[pts[i].best].delta, std::min(pts[i].ub,
                      pts[i].lb[centers[pts[i].best].group_id] - centers[pts[i].best].delta));

                      // the previous best distance becomes 2nd best
                      pts[i].dist_to_center2 = pts[i].ub; 
                      
                      pts[i].best=groups[j].center_ids[k];
                    
                      pts[i].ub = new_d;

                      //log center reassign
                      center_reassignments[i] = 1;
                      //mark centers have changed
                      //yes this is a race, but because we are setting
                      //false to true this is fine
                      centers[pts[i].best].has_changed = true;
                      centers[pts[i].old_best].has_changed=true;

                      
                    }
                    else {
                      if (new_d < pts[i].lb[j]) {
                        pts[i].lb[j] = new_d;
                      }
                    }
                    

                  }
                }
              }
            }
            
            }


      });

      assignment_time = tim.next_time();

      //see how long it takes to get the new_num_members info
      parlay::parallel_for(0,k, [&] (size_t i) {
        centers[i].new_num_members=parlay::filter(pts, [&] (typename ys::point& p) {
          return p.best == i;
        }).size();
      }); // 1 granularity?

      //sanity check
      for (size_t i = 0; i < k; i++) {
        size_t old_mems = parlay::filter(pts, [&] (typename ys::point& p) {
          return p.old_best == i;
        }).size();
        if (old_mems != centers[i].old_num_members) {
          std::cout << "old mems not right " << old_mems << " " << centers[i].old_num_members << std::endl;
          abort();
        }
        if (old_mems != centers[i].new_num_members && centers[i].has_changed==false) {
          std::cout << "center has changed but not marked, aborting " << std::endl;
          abort();
        }
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

#endif //YYSIMP