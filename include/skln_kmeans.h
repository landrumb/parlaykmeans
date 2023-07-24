//Small k large n kmeans algorithm (skin)
//will only run for k < 256
#ifndef SKLN
#define SKLN

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
struct SklnKmeans {

  struct center {
    
    uint8_t id;
    parlay::sequence<float> coordinates;
    parlay::sequence<T> tcoordinates;
    float delta;
    unsigned int old_num_members;
    unsigned int new_num_members;
    //has_changed is true if the center has gained or lost any points
    bool has_changed;

    center(uint8_t id, unsigned short int ds) : id(id) {
      coordinates = parlay::sequence<float>(ds);
      tcoordinates = parlay::sequence<T>(ds);
      delta = std::numeric_limits<float>::max();
      old_num_members = -1;
      new_num_members = -1;
      has_changed = true;
    }


  };


  //compute centers 5

  void compute_centers_compare(T* v , int ns, unsigned short int ds, 
  uint8_t ks, uint8_t* asgs, uint8_t* asgs_old, parlay::sequence<center>& centers, 
  float* center_calc_float, double* center_calc, 
  parlay::sequence<int>& rang) {

    // Compute new centers
    parlay::parallel_for(0,ks,[&](uint8_t j) {
      if (centers[j].has_changed) {
        for (unsigned short int dim = 0; dim < ds; dim++) {
          center_calc[j*ds+dim] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

        }

      }
      
    });


    parlay::sequence<int> changed_points = parlay::filter(rang,[&] (int i) {
      return asgs[i] != asgs_old[i];
    });


    std::cout << "# of reassg: " << changed_points.size() << std::endl;


    parlay::parallel_for(0,ks,[&] (uint8_t j) {
      if (centers[j].has_changed) {
        //TODO wasteful 
        parlay::sequence<int> add_these = parlay::filter(changed_points, [&] (int p) {
        return asgs[p] == j; });
        parlay::sequence<int> sub_these = parlay::filter(changed_points, [&] (int p) {
        return asgs_old[p] == j; });
      //TODO is sequential faster than a reduce? 
      //number of center changes typically very small
      //TODO which order of for loop is better?
      // for (unsigned short int coord = 0; coord < ds; coord++) {
      //   for (size_t ain  = 0; ain < add_these.size(); ain++) {
      //     center_calc[j*ds+coord] += static_cast<double>(v[add_these[ain]*ds+coord]);
      //   }
      //   for (size_t sbin = 0; sbin < sub_these.size(); sbin++) {
      //     center_calc[j*ds+coord] -= static_cast<double>(v[sub_these[sbin]*ds+coord]);
      //   }
        

      // }

      //use a reduce not a sequential for
      for (unsigned short int coord = 0; coord < ds; coord++) {
        center_calc[j*ds+coord] += parlay::reduce(parlay::map(add_these, [&] (int a) {return static_cast<double>(v[a*ds+coord]);})) - 
      parlay::reduce(parlay::map(sub_these, [&] (int a) {return static_cast<double>(v[a*ds+coord]);}));

      }
      


      }
    
      
    }); 

   // t2.next("Added/subbed changes");

    parlay::parallel_for(0,ks,[&](size_t j) {
      
      if (centers[j].new_num_members != 0 && centers[j].has_changed) {
        for (unsigned short int coord = 0; coord < ds; coord++) {
          center_calc[j*ds+coord] /= static_cast<double>(centers[j].new_num_members);
        }
      }
     
    });

    parlay::parallel_for(0,ks,[&] (size_t j) {
      if (centers[j].new_num_members==0 || !centers[j].has_changed) {
        for (unsigned short int coord = 0; coord < ds; coord++) {
          center_calc_float[j*ds+coord] = centers[j].coordinates[coord];
        }
      }
      else {
        for (unsigned short int coord = 0; coord < ds; coord++) {
          center_calc_float[j*ds+coord] = center_calc[j*ds+coord];
        }
      }
    });

    }

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg,
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
    //since we use -1 as a default/empty asg value, k can't be 255
    //in asgs comparing old to new best
    //then again, we could do oldbest = new best-1 for the first iter
    //to avoid this problem
    if (k > 250) {
      std::cout << "k too large, aborting" << std::endl;
      abort();
    }
    uint8_t ks = k;
    if (n >= 1'000'000'001) {
      std::cout << "n too large, aborting" << std::endl;
      abort();
    }
    int ns = n;

    if (d > 1'000) {
      std::cout << "d too big, aborting" << std::endl;
      abort();
    }
    unsigned short int ds = d;

    float setup_time=0;
    float assign_time=0;
    float update_time=0;
    parlay::internal::timer t2;
    t2.start();

    //does using a smaller asgs array actually matter? 
    uint8_t* asgs = new uint8_t[1'000'000'000];
    //_old stores what happened on the previous iteration
    uint8_t* asgs_old = new uint8_t[1'000'000'000];


    //ints safe to use because n <= 1bil, int max is 2bil
    for (int i = 0; i < ns; i++) {
      asgs[i] = asg[i];
    }

    //add initialization stuff here
    parlay::sequence<center> centers = parlay::tabulate(ks,[&] (uint8_t i) {
      return center(i,ds);
    });

    //std::cout << "made it here4" << std::endl;


    for (uint8_t i = 0; i < ks; i++) {
      for (unsigned short int j = 0; j < ds; j++) {
        centers[i].coordinates[j] = *(c + i*ds + j);
        centers[i].tcoordinates[j] = centers[i].coordinates[j];
      }
    }


    //std::cout << "made it here3" << std::endl;


  
    //there must be a better way to get a range
    //TODO ask Laxman
    parlay::sequence<int> rang = parlay::tabulate(ns,[&] (int i) {
      return i;
    });

    //enough room on the stack for these? 
    //not using doubles for now
    float center_calc_float[256000];
    double center_calc[256000];

    //use a histogram to get new num members info, this is faster
    //this is the initial members dist
    auto new_centers_dist = histogram_by_key(parlay::map(rang,[&] (int i) {
      return asg[i];
    }));
    for (uint8_t i = 0; i < ks; i++) {
      centers[new_centers_dist[i].first].new_num_members = new_centers_dist[i].second; //I think this is how it is used? 
      centers[new_centers_dist[i].first].old_num_members = centers[new_centers_dist[i].first].new_num_members;

    }
    

    setup_time = t2.next_time();

          


    std::cout << "made it here5" << std::endl;

    for (size_t iter = 0; iter < max_iter; iter++) {

      std::cout << "made it here7" << std::endl;


      //STEP 1: ASSIGN
      parlay::parallel_for(0, ns, [&](int i) {
        asgs_old[i] = asgs[i]; //set old the number from the previous iter
        auto distances = parlay::delayed::map(centers, [&] (center& q) {
          return D.distance(v+i*ds,q.tcoordinates.begin(),ds);
        });
        asgs[i] = min_element(distances) - distances.begin();
        
      });

      assign_time = t2.next_time();

      std::cout << "made it here 9 " << std::endl;

      //STEP 2: UPDATE

      //STEP 2.1: PREPROCESSING

      //use a histogram to get new num members info, this is faster
      auto new_centers_dist = histogram_by_key(parlay::map(rang,[&] (size_t i) {
        return asgs[i];
      }));


      for (uint8_t i = 0; i < ks; i++) {
        centers[new_centers_dist[i].first].new_num_members = new_centers_dist[i].second; //I think this is how it is used? 
        centers[new_centers_dist[i].first].old_num_members = centers[new_centers_dist[i].first].new_num_members;

      }
      std::cout << "made it here1" << std::endl;

      //update whether or not each center has changed

      //TODO make more efficient? (if possible?) 
      parlay::parallel_for(0,ks,[&] (uint8_t i) {
        centers[i].has_changed = false;
      });

      std::cout << "made it here 40" << std::endl;

      parlay::parallel_for(0,ns,[&] (int i) {
        // if (asgs[i] < 0 || asgs[i] >= ks) {
        //   std::cout << "asgs[" << i << "]=" << asgs[i] << ", aborting" << std::endl;
        //   abort();
        // }
        // if (asgs_old[i] < 0 || asgs_old[i] >= ks) {
        //   std::cout << "asgs_old[" << i << "]=" << asgs_old[i] << ", aborting" << std::endl;
        //   abort();
        // }
        if (asgs[i] != asgs_old[i]) {
          centers[asgs[i]].has_changed = true;
          centers[asgs_old[i]].has_changed = true;
        }
       
      });

      std::cout << "made it here2" << std::endl;

      //STEP 2.2: Actually calculate new centers

      compute_centers_compare(v,ns,ds,ks,asgs,asgs_old,centers,center_calc_float,center_calc,rang);

      //STEP 2.3: Calculate info for next iter, copy over centers

      parlay::sequence<float> deltas = parlay::tabulate(k,[&] (uint8_t i) {
        return D.distance(centers[i].coordinates.begin(),center_calc_float + i*ds,ds);
      });
      float total_diff = *parlay::max_element(deltas); //does max work like this?

      //copy back over centers
      parlay::parallel_for(0,ks,[&](uint8_t i) {
        for (unsigned short int j = 0; j < ds; j++) {
          centers[i].coordinates[j] = center_calc_float[i*ds+j];
          centers[i].tcoordinates[j] = centers[i].coordinates[j];
        }
      });

      update_time=t2.next_time();

      //calculate msse
      //use tcoordinates! (so no casting happens inside)
      //FIXME when tcoordinates not used distances go to 0, 
      //bug in cross type code? 
      double msse = parlay::reduce(parlay::map(rang,[&] (int i) {
        return static_cast<double>(D.distance(v+i*ds,centers[asgs[i]].tcoordinates.begin(),ds));
      }))/ns;

      std::cout << "msse is " << msse << std::endl;


      setup_time += t2.next_time();

      logger.add_iteration(assign_time,update_time,msse,0,0,deltas,setup_time);
      assign_time=0;
      update_time=0;
      setup_time=0;


      if (total_diff <= epsilon) {
        break;
      }


    }
    parlay::parallel_for(0,ns,[&] (int i) {
      asg[i] = asgs[i];
    });
    parlay::parallel_for(0,ks,[&] (uint8_t i) {
      for (unsigned short int coord = 0; coord < ds; coord++) {
        c[i*d+coord] = centers[i].coordinates[coord];
      }
    });


  }
};


#endif //SKLN