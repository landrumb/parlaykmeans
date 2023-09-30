//contains the point, center, and group struct for yy

#ifndef YYSTRUCTS_IMP
#define YYSTRUCTS_IMP

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
struct YyStructsImp {

  //IMP TODO: only an unsigned int needed for the size of data we are looking at
  //might switch to an unsigned int later
  typedef size_t ui; 

  struct point {
    ui best; // the index of the best center for the point
    parlay::slice<T*, T*> coordinates; // the coordinates of the point

    ui id; //an id used for the point, for debugging purposes
    float ub;
    parlay::sequence<float> lb; //lower bounds from a point to each group
    float global_lb; //global lower bound
    ui old_best; //the previous best

    point() : best(-1), coordinates(nullptr, nullptr), 
    ub(std::numeric_limits<float>::max()) {
    }

    point(ui chosen, parlay::slice<T*,T*> coordinates) : best(chosen), 
    coordinates(coordinates.begin(),coordinates.end()), 
    ub(std::numeric_limits<float>::max()) {

    }
  

    point(ui chosen, parlay::slice<T*,T*> coordinates, float ub) : 
    best(chosen), coordinates(coordinates.begin(),coordinates.end()), 
    ub(ub) {

    }

  };

  struct center {
    ui id; // a unique identifier for the center
    ui group_id; //the id of the group that the center belongs to
    parlay::sequence<float> coordinates; // the pointer to coordinates of the center
    float delta;
    ui old_num_members; //how many points belonged to that center
    ui new_num_members; //the number of members this iter
    bool has_changed; //check if the center has changed

    center(ui id, parlay::sequence<float> coordinates) : id(id) {

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
    ui id;
    //store the ids of all the centers belonging to this group
    parlay::sequence<ui> center_ids;

    float max_drift;
  };

};

#endif //YYSTRUCTS_IMP