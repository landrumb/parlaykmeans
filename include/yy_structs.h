//putting the main structs for yy into their own file
//for easy calling by other files

#ifndef YYSTRUCTS
#define YYSTRUCTS

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

template <typename T>
struct YyStructs {


  struct point {

    size_t best; // the index of the best center for the point
    parlay::slice<T*, T*> coordinates; // the coordinates of the point

    size_t id;
    float ub;
    parlay::sequence<float> lb; //lower bounds from a point to each group
    float global_lb; //global lower bound
    float dist_to_center2; //distance to 2nd closest center (exact)
    size_t old_best; //the previous best
    //bool has_changed;

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
    bool has_changed; //check if the center has changed

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

  static void print_point(const point& p) {
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
    

    static void print_group(const group& g) {
      std::cout << "ID: " << g.id << "DEL " << g.max_drift << std::endl;
      for (size_t i = 0; i < g.center_ids.size(); i++) {
        std::cout << g.center_ids[i] << " ";
      }
      std::cout << std::endl;

    }

  static void print_all(parlay::sequence<point> pts, 
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
  static void print41(parlay::sequence<point> pts, 
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
  static void print_target_wide(parlay::sequence<point> pts, 
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
  static void print_target(parlay::sequence<point> pts, 
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


};

#endif //YYSTRUCTS