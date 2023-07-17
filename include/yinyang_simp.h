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
#include "include/utils/union_find.h"

template<typename T>
struct YinyangSimp {

  struct point {

    size_t best; // the index of the best center for the point
    parlay::slice<T*, T*> coordinates; // the coordinates of the point

    size_t id;
    float ub;
    parlay::sequence<float> lb; //lower bounds from a point to each group
    float global_lb; //global lower bound
    float dist_to_center2; //distance to 2nd closest center (exact)
    size_t old_best; //the previous best

    point() : best(-1), coordinates(nullptr, nullptr), 
    ub(std::numeric_limits<float>::max()) {
    }

    point(size_t chosen, parlay::slice<T*,T*> coordinates) : best(chosen), 
    coordinates(coordinates.begin(),coordinates.end()), 
    ub(std::numeric_limits<float>::max()) {

    }

    point(size_t chosen, parlay::slice<T*,T*> coordinates, float ub) : 
    best(chosen), coordinates(coordinates.begin(),coordinates.end()), 
    ub(ub) {

    }

  };

  struct center {
    size_t id; // a unique (hopefully) identifier for the center
    size_t group_id; //the id of the group that the center belongs to
    parlay::sequence<float> coordinates; // the pointer to coordinates of 
    //the center
    float delta;
    size_t old_num_members; //how many points belong to that center
    size_t new_num_members; //the number of members this iter

    center(size_t id, parlay::sequence<float> coordinates) : id(id) {

    this->coordinates = coordinates;
    delta=0;
    group_id = -1;
    old_num_members = 0;
    new_num_members = 0;
    }

    center() : id(-1) {
      delta=0;
      group_id = -1;
      old_num_members = 0;
      new_num_members = 0;

    }

  };
  struct group {
    size_t id;
    //store the ids of all the centers 
    //belonging to this group
    parlay::sequence<size_t> center_ids;

    float max_drift;
  };

  void print_point(const point& p) {
      if constexpr(std::is_same<T,uint8_t>()==true) {
        std::cout << "Po: "<< p.id << std::endl;
        //since the coordinate info for points is fixed, don't need to print
        // for (T* i = p.coordinates.begin(); i != p.coordinates.end(); i = std::next(i) ) {
        //   std::cout << static_cast<int>(*i) << " ";
        // }
        //std::cout << std::endl;
        std::cout << "best: " << p.best << " ";
        std::cout << "ub, 2ub: " << p.ub << " " << p.dist_to_center2 << std::endl;
        std::cout << "global lb; lb: " << p.global_lb << " ";
        for (size_t i = 0; i < p.lb.size(); i++) {
          std::cout << p.lb[i] << " ";

        }
        std::cout << std::endl << std::endl;   

      }
      else {
         std::cout << "Po: "<< std::endl;
        
        for (T* i = p.coordinates.begin(); i != p.coordinates.end(); 
          i = std::next(i) ) {
          std::cout << (*i) << " ";
        }
        std::cout << std::endl;
        std::cout << "best: " << p.best << " ";
        std::cout << "ub, 2ub: " << p.ub << " " << p.dist_to_center2 << std::endl;
        std::cout << "global lb; lb: " << p.global_lb << " ";
        for (size_t i = 0; i < p.lb.size(); i++) {
          std::cout << p.lb[i] << " ";

        }
        std::cout << std::endl; 
      }
      
    }

    void print_center(const center& c) {
      std::cout << "ID" << c.id << std::endl;
      std::cout <<"COORDS ";
      for (size_t i = 0; i < c.coordinates.size(); i++) {
        std::cout << c.coordinates[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Delta: " << c.delta << " Group ID: " << 
      c.group_id << std::endl;

      std::cout << std::endl;
      
    }
    

    void print_group(const group& g) {
      std::cout << "ID: " << g.id << "DEL " << g.max_drift << std::endl;
      for (size_t i = 0; i < g.center_ids.size(); i++) {
        std::cout << g.center_ids[i] << " ";
      }
      std::cout << std::endl;

    }

  void print_all(parlay::sequence<point> pts, 
  parlay::sequence<center> centers, 
  parlay::sequence<group> groups) {
    std::cout << "printing all" << std::endl;
    std::cout << "POINTS: " << std::endl;
    for (size_t i =0; i < pts.size(); i++) {
      print_point(pts[i]);
    }
    std::cout << "CENTERS: " << std::endl;

    for (size_t i = 0; i < centers.size(); i++) {
      print_center(centers[i]);
    }
    std::cout << "GROUPS: " << std::endl;

    for (size_t i = 0; i < groups.size(); i++) {
      print_group(groups[i]);
    }
  }

  //printing all too much junk, so print 41 more useful
  //just pt 41
  void print41(parlay::sequence<point> pts, 
  parlay::sequence<center> centers, 
  parlay::sequence<group> groups, Distance& D) {
    std::cout << "printing 41 all" << std::endl;
    std::cout << "POINTS: " << std::endl;
    print_point(pts[41]);

    std::cout << "dist to center 16: " << 
    pc_dist(pts[41],centers[16],D,pts[41].coordinates.size());
  
    std::cout << "CENTERS: " << std::endl;

    for (size_t i = 0; i < centers.size(); i++) {
      print_center(centers[i]);
    }
    std::cout << "GROUPS: " << std::endl;

    for (size_t i = 0; i < groups.size(); i++) {
      print_group(groups[i]);
    }
  }

  //print a variable pt and all centers this time
  //as well as the distance to the target center
  void print_target_wide(parlay::sequence<point> pts, 
  parlay::sequence<center> centers, 
  parlay::sequence<group> groups, Distance& D, size_t ptarget,
  size_t ctarget) {
    std::cout << "printing 41 all" << std::endl;
    std::cout << "POINTS: " << std::endl;
    print_point(pts[ptarget]);

    std::cout << "dist to center  " << ctarget << ": " <<
    pc_dist(pts[ptarget],centers[ctarget],D,pts[ptarget].coordinates.size())
    << std::endl;
    std::cout << "precision print: " << std::setprecision(10) << 
    pc_dist(pts[ptarget],centers[ctarget],D,pts[ptarget].coordinates.size())
    << " " << pts[ptarget].ub << std::endl;

    std::cout << "CENTERS: " << std::endl;

    for (size_t i = 0; i < centers.size(); i++) {
      print_center(centers[i]);
    }
    std::cout << "GROUPS: " << std::endl;

    for (size_t i = 0; i < groups.size(); i++) {
      print_group(groups[i]);
    }
  }


  //print a variable pt and center this time
  void print_target(parlay::sequence<point> pts, 
  parlay::sequence<center> centers, parlay::sequence<group>& groups, 
  Distance& D, size_t ptarget,
  size_t ctarget) {
    std::cout << "print target call" << std::endl;
    std::cout << "precision 10 " << std::setprecision(10) << std::endl;
    std::cout << "printing target " << ptarget <<  
    " against centers:  " << pts[ptarget].best << ", " << ctarget << std::endl;
    print_point(pts[ptarget]);
    print_center(centers[pts[ptarget].best]);
    print_center(centers[ctarget]);

    std::cout << "dist to center  " << ctarget << ": " <<
    pc_dist(pts[ptarget],centers[ctarget],D,pts[ptarget].coordinates.size())
    << std::endl;

    std::cout << "dist to center  " << pts[ptarget].best << ": " <<
    pc_dist(pts[ptarget],centers[pts[ptarget].best],D,pts[ptarget].coordinates.size())
    << std::endl;
    std::cout << "finished target print>" << 
    std::endl;
   
  }

  //for debugging, get distance between point and center
  float pc_dist(point& p, center& c, Distance& D, size_t d) {
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
  void fill_in_centers(float* c, parlay::sequence<center>& centers, 
  size_t k, size_t d) {
    parlay::parallel_for(0,k, [&] (size_t i) {
    for (size_t j = 0; j < d; j++) {
    centers[i].coordinates[j] = *(c + i*d + j);
    }

    });
  }

  

  //initialize the groups by running a naive kmeans a few times
  void init_groups(float* c, size_t k, size_t d, size_t t, 
    Distance& D, parlay::sequence<center>& centers,
    parlay::sequence<group>& groups) {

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
    group_centers, group_asg,D,logger, 5, 0.0001);
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
  void assert_proper_group_size(parlay::sequence<group>& groups, 
   size_t t, bool DEBUG_FLAG=false) {
    if (DEBUG_FLAG) {
       std::cout << "printing out init groups" << std::endl;
      for (size_t i = 0; i < t; i++) {
        print_group(groups[i]);
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

  //compute centers calculates the new centers
  //puts the new values into the float* center_calc
  //TODO make a different run (no casting) if T is a float
  void compute_centers_filter(
  parlay::sequence<point>& pts, size_t n, size_t d, size_t k, 
  const parlay::sequence<center>& centers, double* center_calc,
  float* center_calc_float) {

    //update centers timer
    parlay::internal::timer t2;
    t2.start();

    std::cout << "using the filter version " << std::endl;

    parlay::sequence<parlay::sequence<size_t>> indices(k);

    t2.next("Made new centers");

    parlay::parallel_for(0,k,[&] (size_t i) {
      auto temp = parlay::filter(pts,[&] (point& p) {
        return p.best==i;
      });
      indices[i] = parlay::map(temp, [&] (point& p) {
        return p.id;
      });
    
    });

    t2.next("Got indices");

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

  

    t2.next("Done with adding");

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    t2.next("Done with copying");
  }


  void compute_centers_comparative(parlay::sequence<point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<center>& centers, 
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


    parlay::sequence<point> changed_points = parlay::filter(pts,[&] (point& p) {
      return p.best != p.old_best;
    });

    t2.next("Filtered");

    parlay::sequence<parlay::sequence<point>> add_these_all(k);
    parlay::sequence<parlay::sequence<point>> sub_these_all(k);
    parlay::parallel_for (0,k,[&] (size_t i) {
      add_these_all[i] = parlay::filter(changed_points, [&] (point& p) {
        return p.best==i;
      });
      sub_these_all[i] = parlay::filter(changed_points, [&] (point& p) {
        return p.old_best == i;
      });

    });

    t2.next("Extra filtering");

    parlay::parallel_for(0,k*d,[&] (size_t jcoord) {
      size_t j = jcoord / d;
      size_t coord = jcoord % d;
      // parlay::sequence<point> add_these = parlay::filter(changed_points, [&] (point& p) {
      //   return p.best==j;
      // });
      // parlay::sequence<point> sub_these = parlay::filter(changed_points, [&] (point& p) {
      //   return p.old_best == j;
      // });
      // for (size_t coord = 0; coord < d; coord++) {
        parlay::sequence<point> add_these = add_these_all[j];
        parlay::sequence<point> sub_these = sub_these_all[j];
        double diff_to_add = 0;
        for (size_t m = 0; m < add_these.size(); m++) {
          diff_to_add += add_these[m].coordinates[coord];
        }
        for (size_t m = 0; m < sub_these.size(); m++) {
          diff_to_add -= sub_these[m].coordinates[coord];
        }
        center_calc[j*d+coord] += diff_to_add;

      // }
      // for (size_t coord = 0; coord < d; coord++) {
      //   center_calc[j*d+coord] = center_calc[j*d+coord] + 
      //   parlay::reduce(parlay::map(add_these,[&] (point& p) {return static_cast<double>(p.coordinates[coord]);} ),100) - 
      //   parlay::reduce(parlay::map(sub_these, [&] (point& p) {return static_cast<double>(p.coordinates[coord]);}),100);

      // }
    });


    // //Remark: changed_points usually very small so this is okay
    // for (size_t i = 0; i < changed_points.size(); i++) {
    //   //TODO should this inner loop be ||? (on d)
    //   point p = pts[i];
    //   for (size_t coord =0 ; coord < d; coord++) {
      //TODO CAST TO DOUBLE
    //     center_calc[p.best * d + coord] += p.coordinates[coord];
    //     center_calc[p.old_best *d + coord] -= p.coordinates[coord];
    //   }

    // }


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

  //computes the drift, then initializes the new centers
  float update_centers_drift(parlay::sequence<point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<center>& centers, Distance& D, 
  parlay::sequence<group>& groups, size_t t, float* center_calc_float) {
   
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

  void set_point_global_lb(point& p, parlay::sequence<group>& groups,
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

  void debug_group_asg_centers(parlay::sequence<center>& centers,float* group_centers, float* group_asg,
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

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {

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
    parlay::sequence<point> pts = parlay::tabulate<point>(n, [&] (size_t i) {
    return point(asg[i],parlay::slice(v+i*d, v+i*d + d));

    });

    parlay::parallel_for(0,n,[&] (size_t i) {
      pts[i].id=  i;
    });

    //create the centers
    parlay::sequence<center> centers = parlay::tabulate<center>(k, [&] (size_t i) {
    return center(i, parlay::sequence<float>(d));
    });

    //fill in the centers
    fill_in_centers(c,centers,k,d);
    
    //making sure t is at least 1
    //TODO change back to k/10 (k/20 loads faster)
    size_t t = std::max((size_t) k/20,(size_t) 1); 
    
  
    parlay::sequence<group> groups(t);
    //initialize the groups
    init_groups(c,k,d,t,D,centers,groups);

    logger.add_iteration(0,0,45,0,0,
    parlay::sequence<float>(k,0),tim.next_time());

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
        
        auto distances = parlay::delayed::map(centers, [&](center& q) {
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

    logger.add_iteration(0,0,48,0,0,
    parlay::sequence<float>(k,0),tim.next_time());

    //set num_members for the centers
    parlay::parallel_for(0,k,[&] (size_t i) {
      auto belonging_pts = parlay::filter(pts,[&] (point& p) {
        return p.best==i;
      });
      centers[i].new_num_members=belonging_pts.size();
      centers[i].old_num_members=belonging_pts.size();
    });

    logger.add_iteration(0,0,49,0,0,parlay::sequence<float>(k,0),
    tim.next_time());

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
    //Step 3: Repeat until convergence

    //for center calculation
    double* center_calc = new double[k*d];
    float* center_calc_float = new float[k*d];

    bool first_time = true;

    
    while (true) {
      //the first iteration, we don't have previous centers to compare to
      //so we must do a full-adding version of compute centers
      if (first_time) {
        compute_centers_filter(pts,n,d,k,centers,center_calc,center_calc_float);
        first_time=false;
      }
      else {
        compute_centers_filter(pts,n,d,k,centers,center_calc,center_calc_float);
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
      logger.add_iteration(assignment_time,update_time,squared_error,parlay::reduce(distance_calculations),
      parlay::reduce(center_reassignments),deltas);

      //convergence check
      if (iters >= max_iter || total_diff <= epsilon) break;

      iters += 1;
      distance_calculations = parlay::sequence<size_t>(n,0); //reset distance calc counter
      center_reassignments = parlay::sequence<uint8_t>(n,0);

      std::cout << "iter: " << iters << std::endl;

      //3.2: Group filtering

      parlay::parallel_for(0,n,[&](size_t i) {
        //update bounds
        pts[i].ub += centers[pts[i].best].delta; 
        
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
                      centers[groups[j].center_ids[k]].new_num_members += 1;
                      centers[pts[i].best].new_num_members -= 1;

                      pts[i].old_best = pts[i].best; //save the previous best in old_best

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