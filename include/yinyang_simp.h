//yinyang but no local filtering
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

template<typename T>
struct YinyangSimp {

  bool DEBUG_VD = false;
  size_t TARGET_POINT = 1971;
  size_t TARGET_CENTER = 17;

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

  //print a variable pt and center this time
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

  //compute centers calculates the new centers
  //returns: a sequence of centers
  //TODO make a different run (no casting) if T is a float
  parlay::sequence<center> compute_centers_vd(
  const parlay::sequence<point>& pts, size_t n, size_t d, size_t k, 
  const parlay::sequence<center>& centers) {

    parlay::sequence<parlay::sequence<size_t>> indices(k);

    parlay::sequence<center> new_centers(k);
    for (size_t i = 0; i < k; i++) {
    new_centers[i].id = i;
    new_centers[i].coordinates=parlay::sequence<float>(d,0);
    }


    for (size_t i = 0; i < n; i++) {

    indices[pts[i].best].push_back(i); //it's called best not id!!!

    }
    if (DEBUG_VD) {
    std::cout << "Debugging: printing out center counts:\n";
    for (size_t i = 0; i < k; i++) {
    std::cout << indices[i].size() << std::endl;
    }

    }


    parlay::parallel_for (0, k*d, [&] (size_t icoord){
    size_t i = icoord / d;
    size_t coord = icoord % d;


    //if there are no values in a certain center, just keep the center 
    //where it is
    if (indices[i].size() > 0) { //anti_overflow_avg or reduce?? 
    //note the static cast to double here, because points are whatever
    new_centers[i].coordinates[coord] = static_cast<float>(reduce(parlay::map(indices[i],[&] 
    (size_t ind) {return static_cast<double>(
    pts[ind].coordinates[coord]);})) / indices[i].size()); 
    //normal averaging now

    }
    else { 
    new_centers[i].coordinates[coord] = centers[i].coordinates[coord];
    }


    });


    return new_centers;


  }

  //given a point and a list of distances from the point to the closest
  //and second closest center in each group, set p's initial
  //best center, ub, and distance to 2nd closest center

  //at index i:
  //<size_t id, float dist>
  //id: id of center in group i that is closest to p
  //dist: distance between center id and p 
  void point_set_best_ub(point& p, 
  parlay::sequence<std::pair<size_t,float>>& distances,
  parlay::sequence<std::pair<size_t,float>>& distances2nd,
  parlay::sequence<center>& centers) {
    size_t best_index = -1;
    float best_dist = std::numeric_limits<float>::max();
    float second_dist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < distances.size(); i++) {
      if (distances[i].second < best_dist) {
        best_index = i;
        second_dist = best_dist;
        best_dist = distances[i].second;
      }
      else if (distances[i].second < second_dist) {
        second_dist = distances[i].second;
      }
      //no need to compare distances2nd against the best dist
      if (distances2nd[i].second < second_dist) {
        second_dist = distances2nd[i].second;
      }

    }
    p.best = distances[best_index].first;
    p.old_best = p.best;
    p.ub = distances[best_index].second;
    p.dist_to_center2 = second_dist;
    centers[p.best].old_num_members += 1; //keep track of membership

  }

  //TODO: If the type is float, then run a separate branch, 
  //no need to convert to float

  //find the distance from a point to each group (closest and 2nd closest)
  //records a pair <size_t, float>: the id of the center giving a 
  //certain distance, and the distance itself
  void distances_to_groups(const point&p, parlay::sequence<center>& centers,
  parlay::sequence<group>& groups,
  Distance& D, size_t t, size_t d,
  parlay::sequence<std::pair<size_t,float>>& distances, 
  parlay::sequence<std::pair<size_t,float>>& distances2nd) {

    //convert point coordinates to float buffer
    float buf[2048];
    T* it = p.coordinates.begin();
    for (size_t i = 0; i < d; i++) {
    buf[i]=*it;
    it += 1; //add 1 for next?
    }

    //for each group
    parlay::parallel_for(0,t,[&] (size_t i) {
    //calculate the distance from each group member to the point
    auto dist = parlay::map(groups[i].center_ids, 
    [&](size_t j) { //note j is from centers because it's a map not a for loop (error? TODO)
    return std::make_pair(j,
    sqrt_dist(buf, make_slice(centers[j].coordinates).begin(),d,D));
    });

    //if a group only has only member, then 2nd should store a very large distance
    if (groups[i].center_ids.size() == 1) {
      distances[i] = dist[0];
      distances2nd[i] = std::make_pair(0,std::numeric_limits<float>::max());

    }
    else {
      //find the closest and 2nd closest center
      //TODO make parallel
      size_t min_index = -1;
      size_t min2_index = -1; //2nd smallest
      if (dist[0].second < dist[1].second) {
      min_index = 0;
      min2_index = 1;
      }
      else {
      min_index = 1;
      min2_index = 0;
      }
      for (size_t j = 2; j < groups[i].center_ids.size(); j++) {
      if (dist[j].second < dist[min2_index].second) {
      if (dist[j].second < dist[min_index].second) {
      min2_index = min_index;
      min_index = j;

      }
      else {
      min2_index = j;
      }
      }

      }
      distances[i] = dist[min_index];
      distances2nd[i] = dist[min2_index];
      

    }
    });





  }

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
  float* group_centers, size_t* group_asg, Distance& D, 
  parlay::sequence<center>& centers) {

    std::cout << "Debugging init groups " << std::endl;
   
    LazyStart<float> init;
    init(c,k,d,t,group_centers,group_asg,D);
    NaiveKmeans<float> run;
    run.cluster(c,k,d,t,
    group_centers, group_asg,D, 5, 0.0001);


    parlay::parallel_for(0,k,[&] (size_t i) {
      centers[i].group_id = group_asg[i];
    });


  }

  //note that groups can have different sizes!
  //now the groups can have any positive size, even size 1
  void assert_proper_group_size(parlay::sequence<group>& groups, 
   size_t t) {
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

  //on the first iteration, for all points, set the point's best center
  //and the upper and lower bounds, and the 2nd closest center
  //(all the init point info)
  void init_all_point_bounds(parlay::sequence<point>& pts, size_t n, 
  parlay::sequence<parlay::sequence<std::pair<size_t,float>>>& distances, 
  parlay::sequence<parlay::sequence<std::pair<size_t,float>>>& distances2nd,
  parlay::sequence<center>& centers, size_t t) {
    parlay::parallel_for(0,n,[&] (size_t i) {

      //assign closest and 2nd closest center
      point_set_best_ub(pts[i],distances[i],distances2nd[i], centers);

    // pts[i].best=min_element(distances[i])-distances[i].begin();
    // pts[i].ub = distances[pts[i].best];
    //could do in parallel? TODO
    for (size_t j = 0; j < t; j++) {
      //in the group containing the best, you take the 2nd best from
      //the group
      if (centers[pts[i].best].group_id==j) {
        pts[i].lb[j] = distances2nd[i][j].second;
      }
      //but in any other group, the lower bound is the group element
      //closest to p
      else {
        pts[i].lb[j] = distances[i][j].second;
      }

    }


    });

  }

  float update_centers_drift(parlay::sequence<point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<center>& centers, Distance& D, 
  parlay::sequence<group>& groups, size_t t) {
    // Compute new centers
    parlay::sequence<center> new_centers = compute_centers_vd(pts, n, d, k, 
    centers);

    // Check convergence
    float total_diff = 0.0;
    for (size_t i = 0; i < k; i++) { 
      new_centers[i].delta = sqrt_dist(centers[i].coordinates.begin(), 
      new_centers[i].coordinates.begin(),d,D); // lol not a k here, a d
      new_centers[i].group_id = centers[i].group_id; //need to copy over this info!
      total_diff += new_centers[i].delta;
    }

    centers = std::move(new_centers); //Q* cause of error? TODO move failing perhaps no evidence though

    //for each group, get max drift for group
    parlay::parallel_for(0,t,[&] (size_t i) {
      auto drifts = parlay::map(groups[i].center_ids, [&] (size_t j) {
      return centers[j].delta;

      });
      
      groups[i].max_drift = *max_element(drifts);

    });

    return total_diff;

  }

  float update_centers_drift_comparative(parlay::sequence<point>& pts, size_t n, size_t d, 
  size_t k, parlay::sequence<center>& centers, Distance& D, 
  parlay::sequence<group>& groups, size_t t) {
    // Compute new centers
    double* center_calc = new double[k*d]; //hold new centers here for a bit
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      center_calc[i] = centers[j].coordinates[dim] * centers[j].old_num_members;

    });
    for (size_t i = 0; i < n; i++) {
      if (pts[i].best != pts[i].old_best) {
        for (size_t j = 0; j < d; j++) {
          center_calc[pts[i].best * d + j] += pts[i].coordinates[j];
          center_calc[pts[i].old_best *d + j] -= pts[i].coordinates[j];
        }

      }
      
    }
    parlay::parallel_for(0,k*d,[&](size_t i) {
      size_t j = i / d;
      center_calc[i] /= centers[j].new_num_members;
    });

    float* center_calc_float = new float[k*d];
    parlay::parallel_for(0,k*d,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });
  
   

    // Check convergence
    float total_diff = 0.0;
    for (size_t i = 0; i < k; i++) { 
      centers[i].delta = sqrt_dist(centers[i].coordinates.begin(), 
      center_calc_float + i*d,d,D); // lol not a k here, a d
      total_diff += centers[i].delta;
    }

    parlay::parallel_for(0,k*d,[&] (size_t i) {
      size_t j = i / d;
      size_t dim = i % d;
      centers[i].coordinates[j] = center_calc[i];
      
    });

    //for each group, get max drift for group
    parlay::parallel_for(0,t,[&] (size_t i) {
      auto drifts = parlay::map(groups[i].center_ids, [&] (size_t j) {
      return centers[j].delta;

      });
      
      groups[i].max_drift = *max_element(drifts);

    });

    delete[] center_calc;
    delete[] center_calc_float;

    return total_diff;

  }

  void set_point_global_lb(point& p, parlay::sequence<group>& groups,
  size_t t) {
    p.global_lb = std::numeric_limits<float>::max();
    for (size_t j = 0; j < t; j++) {
      //max not working? so using comparison
      //not working p.lb[j] = std::max(0.0,p.lb[j]-groups[j].max_drift);
      if (p.lb[j] - groups[j].max_drift > 0) {
        p.lb[j] = p.lb[j] - groups[j].max_drift;
      }
      else {
        p.lb[j] = 0;
      }

      //reduce the global lower bound if possible
      if (p.global_lb > p.lb[j]) {
        p.global_lb = p.lb[j];
      } 
    }

  }

  //yinyang needs sqrt dist for triangle inequality to work
  float sqrt_dist(float* a, float* b, size_t d, Distance& D) {
    return std::sqrt(D.distance(a,b,d));
  }


  

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
  Distance& D, size_t max_iter, double epsilon) {


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


    //format the data according to our naive run
    parlay::sequence<point> pts = parlay::tabulate<point>(n, [&] (size_t i) {
    return point(asg[i],parlay::slice(v+i*d, v+i*d + d));

    });
    //TODO remove pt.id purely for debugging
    parlay::parallel_for(0,n,[&] (size_t i) {
      pts[i].id=  i;
    });

    //create the centers
    parlay::sequence<center> centers = parlay::tabulate<center>(k, [&] (size_t i) {
    return center(i, parlay::sequence<float>(d));
    });

    //fill in the centers
    fill_in_centers(c,centers,k,d);

    // std::cout << "printing out init points, centers";
    // for (size_t i = 0; i < n; i++) {
    //   print_point(pts[i]);
    // }
    // for (size_t i = 0; i < k; i++) {
    //   print_center(centers[i]);
    // }
    
    
    size_t t = k/10; //TODO change back to k/10
    
    //cluster on the groups initially using NaiveKmeans
    float* group_centers = new float[t * d];
    size_t* group_asg = new size_t[k];

    //initialize the groups
    init_groups(c,k,d,t,group_centers,group_asg,D,centers);
    //these checks no longer needed
    // std::cout << "Couting same info as above (group asg) but outside of the function" << std::endl;
    // std::cout << "Group asg: " << std::endl;
    // for (size_t i = 0; i < k; i++) {
    //   std::cout << group_asg[i] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Couting group centers: " << std::endl;
    // for (size_t i = 0; i < t; i++) {
    //   for (size_t j = 0; j < d; j++) {
    //     std::cout << group_centers[i*d+j] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << "Couting group ids" << std::endl;
    // for (size_t i = 0; i < k; i++) {
    //   std::cout << centers[i].group_id << " ";
    // }
    // std::cout << std::endl;
    

    parlay::sequence<group> groups(t);

    //need to set group's id!! lol
    for (size_t i = 0; i < t; i++) {
      groups[i].id = i;
    }

    


    //Sadly 
    //sequential assigning the groups their centers
    for (size_t i = 0; i < k; i++) {
      groups[centers[i].group_id].center_ids.push_back(i);
      
    }

     std::cout << "printing out init groups" << std::endl;
    for (size_t i = 0; i < t; i++) {
      print_group(groups[i]);
    }

    //confirm group assignment happened properly
    //(checking a necessary not sufficient condition)
    assert_proper_group_size(groups,t);

    //Yinyang first does a normal iteration of kmeans:   

    //Assignment
    //distance from each group to a point
    //at index (i,j) is the the id of the closest center in group j,
    //as well as the distance between point i and the closest center
    parlay::sequence<parlay::sequence<std::pair<size_t,float>>> 
      distances(n,parlay::sequence<std::pair<size_t,float>>(t));
    //distance from 2nd closest point in each group to a point
    parlay::sequence<parlay::sequence<std::pair<size_t,float>>> 
      distances2nd(n,parlay::sequence<std::pair<size_t,float>>(t));
    
    // Assign each point to the closest center
    parlay::parallel_for(0, pts.size(), [&](size_t i) {
    distances_to_groups(pts[i],centers,groups,D,t,d,distances[i],
    distances2nd[i]);


    pts[i].lb = parlay::sequence<float>(t); //initialize the lower bound sequence
    });


    //print out distances, distance2nd

    std::cout << "printing out distances, distances2nd for pt 0 " << std::endl;
    for (size_t i = 0; i < t; i++) {
      std::cout << distances[0][i].first << " " << distances[0][i].second << 
      " | " << distances2nd[0][i].first << " " << distances2nd[0][i].second << std::endl;
    }

    //actually take best
    init_all_point_bounds(pts, n, distances, distances2nd, centers, t);

    parlay::parallel_for(0,k,[&] (size_t i) {
      centers[i].new_num_members = centers[i].old_num_members;
    });

    //debugging
    // {
    // size_t elt_counter = 0;
    // for (size_t i = 0; i < n; i++) {
    //   elt_counter += centers[i].old_num_members;

    // }
    // if (elt_counter != n) {
    //   std::cout << "error in num_members assignment: " << elt_counter << std::endl;
    //   abort();
    // }}

    //iters must start at 1!! (was like that a couple iterations ago
    //ironically)
    size_t iters = 1; 
    float total_diff = 0.0;

    //confirm that the initial point assigning went correctly
    //print_all(pts,centers,groups);
    //experiment: are we setting the lbs initially wrong, or 
    //are they becoming wrong over time?
    // for (int i = 0; i < n; i++) {
    //   pts[i].global_lb = -1;
    //   for (size_t j = 0; j < t; j++) {
    //     pts[i].lb[j] = -1;
    //   }

    // }


    //Step 3: Repeat until convergence

    
    while (true) {
      // std::cout << "print center/point before drift" << std::endl;
      // print_target(pts,centers,groups,D,TARGET_POINT,TARGET_CENTER);


      
      //3.1: update centers
      //TODO use the yinyang fast update centers method (currently doing
      //a naive update centers)
      //comparative is the fast method
      total_diff = update_centers_drift(pts,n,d,k,centers,D,groups, t);

      //convergence check
      if (iters >= max_iter || total_diff < epsilon) break;

      iters += 1;

      
      // std::cout << "print 41 after drift" << std::endl;
      // print_target(pts,centers,groups,D,TARGET_POINT,TARGET_CENTER);

       
      

      //check everything after each iter
      std::cout << "iter: " << iters << std::endl;

      //3.2: Group filtering

      //update bounds
      parlay::parallel_for(0,n,[&](size_t i) {
        pts[i].ub += centers[pts[i].best].delta; 
        
        set_point_global_lb(pts[i],groups,t);
        // if (i==TARGET_POINT) {
        //   std::cout << "printing point 41 after set global lb iter " 
        // << iters << std::endl;
        //     print_target(pts,centers,groups,D,TARGET_POINT,TARGET_CENTER);


        // }
        
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
          //do local filtering (3.3)
          //TODO add in local filtering 
          //doing the brute force calculations for now

          //for each group
          //!! in this for loop we are using size_t * n memory
          //anyway (wasteful?) (oh not *d so still okay)
          //overall less memory than our points take
          for (size_t j = 0; j < t; j++) {
            //if the lower bound of group j is smaller than the 
            //upper bound of pt i
            if (pts[i].ub > pts[i].lb[j]) {
              //look at each group member (ie each center)

              //reset the lower bound, make it smallest distance we calculate
              //that's not the real distance away
              pts[i].lb[j] = std::numeric_limits<float>::max(); 

              for (size_t k = 0; k < groups[j].center_ids.size(); k++) {

                
                //since we skip the local filter, we do a distance calculation
                float new_d = sqrt_dist(
                buf,
                parlay::make_slice(
                  centers[groups[j].center_ids[k]].coordinates).begin(),
                d,D);

                //note that the ub is tight rn,
                //that ub IS the distance to the previously 
                //closest center
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
                    pts[i].ub - centers[pts[i].best].delta;

                  // the previous best distance becomes 2nd best
                  pts[i].dist_to_center2 = pts[i].ub; 
                  
                  //and we can also use this 2nd best distance as a lower
                  //bound for group j, as the 2nd best distance is 
                  //necessarily better than all of the other elements
                  //in group j (even though elt k is THE best)
                  //without local filtering, don't need this
                  // if (pts[i].ub < pts[i] < lb[j]) {
                  //   pts[i].lb[j] = pts[i].ub;

                  // }
                    
                  pts[i].best=groups[j].center_ids[k];
                
                  pts[i].ub = new_d;
                  
                  
                }
                else {
                  if (new_d < pts[i].lb[j]) {
                    pts[i].lb[j] = new_d;
                  }
                }
                

              }
            }
            // else { //unsure why this is needed? TODO try removing
            //   pts[i].lb[j] -= groups[j].max_drift;
            // }
            

          }
        }
        
      
        }


      });



    }
    // std::cout << "post while loop printing, looking at bounds " << std::endl;
    // print_target(pts,centers,groups,D,TARGET_POINT,TARGET_CENTER);


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


    delete[] group_centers; //memory cleanup
    delete[] group_asg;

  }

};

#endif //YYSIMP