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
struct NiskKmeans {

        struct center {
    
    uint8_t id;
    parlay::sequence<float> coordinates;
    parlay::sequence<T> tcoordinates;
    float delta;
    unsigned int old_num_members;
    unsigned int new_num_members;
    //has_changed is true if the center has gained or lost any points
    bool has_changed;
    bool stable;
    float radius;
    parlay::sequence<uint8_t> neighbors;
    parlay::sequence<uint8_t> old_neighbors;
    
    center(uint8_t id, unsigned short int ds) : id(id) {
      coordinates = parlay::sequence<float>(ds);
      tcoordinates = parlay::sequence<T>(ds);
      delta = std::numeric_limits<float>::max();
      old_num_members = -1;
      new_num_members = -1;
      has_changed = true;
      stable = false;
    }


  };

  //the first iter needs its own compute center because
  //we can't compare if the assignments are random like Lazy init does
  //compute centers calculates the new centers
  //puts the new values into the float* center_calc
  //TODO make a different run (no casting) if T is a float
  //typename needed because point a dependent type
  static void compute_centers_filter(
  T* v, int ns, unsigned short int ds, uint8_t ks, uint8_t* asgs, parlay::sequence<center>& centers, float* center_calc_float, double* center_calc, parlay::sequence<int>& rang) {

    //update centers timer
    //commenting out timer prints because too much printing
    parlay::internal::timer t2;
    t2.start();

    //TODO pass this in to save time?
     auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rang,[&] (int i) {
        return std::pair(asgs[i],i);
      }));


    //t2.next("Got indices");

    parlay::parallel_for (0, ks, [&] (size_t j){
      size_t picked_center = pts_grouped_by_center[j].first;
      
   

    //if there are no values in a certain center, just keep the center 
    //where it is
    //need to cast to double inside, before the reduce starts adding
    if (pts_grouped_by_center[j].second.size() > 0) { 
      for (size_t coord = 0; coord < ds; coord++) {
        center_calc[picked_center*ds+coord] = static_cast<double>(parlay::reduce(parlay::map(pts_grouped_by_center[j].second, [&]
        (size_t ind) {
          return static_cast<double>(v[ind*ds+coord]);
        })))/pts_grouped_by_center[j].second.size();
      }
  
    }
    else { 
      for (size_t coord = 0; coord < ds; coord++) {
        center_calc[picked_center*ds + coord] = centers[picked_center].coordinates[coord];

      }
    }

    }); 

  

    //t2.next("Done with adding");

    parlay::parallel_for(0,ks*ds,[&] (size_t i) {
      center_calc_float[i] = center_calc[i];
    });

    //t2.next("Done with copying");
  }


  //compute centers 5

  void compute_centers_compare(T* v , int ns, unsigned short int ds, 
  uint8_t ks, uint8_t* asgs, uint8_t* asgs_old, parlay::sequence<center>& centers, 
  float* center_calc_float, double* center_calc, 
  parlay::sequence<int>& rang) {

    // Compute new centers
    parlay::parallel_for(0,ks,[&](size_t j) {
      if (centers[j].has_changed) {
        for (unsigned short int dim = 0; dim < ds; dim++) {
          center_calc[j*ds+dim] = static_cast<double>(centers[j].coordinates[dim]) * centers[j].old_num_members;

        }

      }
      
    });


    parlay::sequence<int> changed_points = parlay::filter(rang,[&] (int i) {
      return asgs[i] != asgs_old[i];
    });


    //std::cout << "# of reassg: " << changed_points.size() << std::endl;


    parlay::parallel_for(0,ks,[&] (uint8_t j) {
      if (centers[j].has_changed) {
        //TODO wasteful 
        parlay::sequence<int> add_these = parlay::filter(changed_points, [&] (int p) {
        return asgs[p] == j; });
        parlay::sequence<int> sub_these = parlay::filter(changed_points, [&] (int p) {
        return asgs_old[p] == j; });

      //use a reduce not a sequential for
      for (unsigned short int coord = 0; coord < ds; coord++) {
        center_calc[j*ds+coord] = center_calc[j*ds+coord] + parlay::reduce(parlay::map(add_these, [&] (int a) {return static_cast<double>(v[a*ds+coord]);})) - 
      parlay::reduce(parlay::map(sub_these, [&] (int a) {return static_cast<double>(v[a*ds+coord]);}));

      }
      


      }
    
      
    }); 

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
    //sort neighbor id
    void sort_neighbor(parlay::sequence<center>& centers, 
    parlay::sequence<parlay::sequence<float>>& dist_matrix, uint8_t id) {

      parlay::sequence<uint8_t>& neighbors = centers[id].neighbors;
      //std::cout << "# of neighbors of center " << id << ": " << centers[id].neighbors.size() << std::endl;
      int min_index;
      float min_val;
      if(neighbors.size() == 0){
        return;
      }
      for (size_t i = 0 ; i < neighbors.size()-1; i++) {
        //std::cout << "i: " << static_cast<int>(i) << std::endl; 
        min_index = i;
        min_val = dist_matrix[id][neighbors[i]];
        for (size_t j = i+1; j < neighbors.size(); j++) {
          if (dist_matrix[id][neighbors[j]] < min_val) {
            min_val = dist_matrix[id][neighbors[j]];
            min_index = j;
          }
          
        }

        uint8_t temp = neighbors[i];
        neighbors[i] = neighbors[min_index];
        neighbors[min_index] = temp;
      }
     //std::cout << "done with sort neighbor" << std::endl;
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

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    /////Part 1: SETUP
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

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
    //plug the assignments into asgs, which is what we will be using
    for (int i = 0; i < ns; i++) {
      asgs[i] = asg[i];
    }

    //add initialization stuff here
    //initialize the centers
    parlay::sequence<center> centers = parlay::tabulate(ks,[&] (uint8_t i) {
      return center(i,ds);
    });

    //std::cout << "made it here4" << std::endl;

    //fill in the centers with the coordinates given in c
    for (size_t i = 0; i < ks; i++) {
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
    float center_calc_float[256000];
    double center_calc[256000];

    //use a histogram to get new num members info, this is faster
    //this is the initial members dist
    auto new_centers_dist0 = histogram_by_key(parlay::map(rang,[&] (int i) {
      return asgs[i];
    }));

    if (new_centers_dist0.size() != ks) {
      std::cout << "center lost all points, aborting" << std::endl;
      abort();
    }

    for (uint8_t i = 0; i < ks; i++) {
      centers[new_centers_dist0[i].first].new_num_members = new_centers_dist0[i].second; //I think this is how it is used? 
      centers[new_centers_dist0[i].first].old_num_members = centers[new_centers_dist0[i].first].new_num_members;

    }
    //dist_matrix is a 2D list that gives the distance from center i to center j at index (i,j)
    parlay::sequence<parlay::sequence<float>> dist_matrix(ks,parlay::sequence(ks,std::numeric_limits<float>::max()));

    setup_time = t2.next_time();

        
    //std::cout << "made it here5" << std::endl;

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    /////Part 2: ITERATIONS: The zeroth iter
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
  
    // std::cout << "printing first 10 coords of each center: " << std::endl;
    //   for (uint8_t i = 0; i < ks; i++) {
    //     for (int j = 0; j < 10; j++) {
    //       std::cout << centers[i].coordinates[j] << " ";
    //     }
    //     std::cout << std::endl;
    //   }
    //STEP 1: ASSIGN
    parlay::parallel_for(0, ns, [&](int i) {
      asgs_old[i] = asgs[i]; //set old the number from the previous iter
      auto distances = parlay::delayed::map(centers, [&] (center& q) {
        return D.distance(v+i*ds,q.tcoordinates.begin(),ds);
      });
      asgs[i] = min_element(distances) - distances.begin();
      
    });
      
    assign_time = t2.next_time();

    //STEP 2: UPDATE

    //STEP 2.1: PREPROCESSING


    //use a histogram to get new num members info, this is faster
    auto new_centers_dist = histogram_by_key(parlay::map(rang,[&] (int i) {
      return asgs[i];
    }));

    //instead of catching error down the line catch it here
    if (new_centers_dist.size() != ks) {
      std::cout << "center not represented aborting" << std::endl;
      abort();
    }

    //update centers new_num_members and old_num_members
    for (uint8_t i = 0; i < ks; i++) {

 
      centers[new_centers_dist[i].first].new_num_members = new_centers_dist[i].second; //I think this is how it is used? 
      //centers[new_centers_dist[i].first].old_num_members = centers[new_centers_dist[i].first].new_num_members;

    }
    //std::cout << "made it here1" << std::endl;

    //update whether or not each center has changed
    //TODO make more efficient? (if possible?) 
    parlay::parallel_for(0,ks,[&] (uint8_t i) {
      centers[i].has_changed = false;
    });

    //std::cout << "made it here 40" << std::endl;

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

    //std::cout << "made it here2" << std::endl;

    //STEP 2.2: Actually calculate new centers

    compute_centers_filter(v,ns,ds,ks,asgs,centers,center_calc_float,center_calc,rang);

    //STEP 2.3: Calculate info for next iter, copy over centers
    //calculate how much each center moved
    parlay::sequence<float> deltas = parlay::tabulate(k,[&] (size_t i) {
        centers[i].delta = D.distance(centers[i].coordinates.begin(),center_calc_float + i*ds,ds);
        return centers[i].delta;
      });
    float total_diff = *parlay::max_element(deltas); //does max work like this?

    //copy back over centers
    parlay::parallel_for(0,ks,[&](size_t i) {
      for (unsigned short int j = 0; j < ds; j++) {
        centers[i].coordinates[j] = center_calc_float[i*ds+j];
        centers[i].tcoordinates[j] = centers[i].coordinates[j];
      }
      centers[i].old_num_members = centers[i].new_num_members; //do need this set afterwards
    });

    update_time=t2.next_time();


    //calculate msse
    //use tcoordinates! (so no casting happens inside)
    //FIXME when tcoordinates not used distances go to 0, 
    //bug in cross type code? 
    double msse = parlay::reduce(parlay::map(rang,[&] (int i) {
      return static_cast<double>(D.distance(v+i*ds,centers[asgs[i]].tcoordinates.begin(),ds));
    }))/ns;

    //std::cout << "msse is " << msse << std::endl;

    setup_time += t2.next_time();

    logger.add_iteration(assign_time,update_time,msse,0,0,deltas,setup_time);
    assign_time=0;
    update_time=0;
    setup_time=0;


    if (total_diff <= epsilon) {
      std::cout << "finishing after one iter" << std::endl;
      //copy back data points
      parlay::parallel_for(0,ns,[&] (int i) {
        asg[i] = asgs[i];
      });
      parlay::parallel_for(0,ks,[&] (size_t i) {
        for (unsigned short int coord = 0; coord < ds; coord++) {
          c[i*d+coord] = centers[i].coordinates[coord];
        }
      });
      return; //leave
    }
    //DEBUGGING
    // std::cout << "after iter 0: " << std::endl;
    // int debug_num_members = 0;
    // int debug_num_members_new = 0;
    // for (uint8_t i = 0; i < ks; i++) {
    //   std::cout << "center " << static_cast<int>(i) << "has " << centers[i].old_num_members << std::endl;
    //   debug_num_members += centers[i].old_num_members;
    //   debug_num_members_new += centers[i].new_num_members;
    // }
    // if (debug_num_members != ns || debug_num_members_new != ns) {
    //   std::cout << "members don't add up, aborting" << std::endl;
    //   std::cout << debug_num_members << " " << debug_num_members_new << std::endl;
    //   abort();
    // }

  
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    /////Part 3: ITERATIONS (All of the other iterations)
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (size_t iter = 1; iter < max_iter; iter++) {

      std::cout << "Starting iter " << iter << std::endl;

      // std::cout << "printing first 10 coords of each center: " << std::endl;
      // for (uint8_t i = 0; i < ks; i++) {
      //   for (int j = 0; j < 10; j++) {
      //     std::cout << centers[i].coordinates[j] << " ";
      //   }
      //   std::cout << std::endl;
      // }
      // std::cout << "printing first 10 tcoords of each center: " << std::endl;
      // for (uint8_t i = 0; i < ks; i++) {
      //   for (int j = 0; j < 10; j++) {
      //     std::cout << static_cast<int>(centers[i].tcoordinates[j]) << " ";
      //   }
      //   std::cout << std::endl;
      // }

      //DEBUG TEST ASGS
      // std::cout << "printing 100 asgs:" << std::endl;
      // for(int i = 0; i < 100; i++){
      //   std::cout << static_cast<int>(asgs[i]) << std::endl;
      // }

      //update radii
      //holds pairs of center ids and the points associated with that center
      auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rang,[&] (int i) {
        return std::pair(asgs[i],i);
      }));

      //std::cout << "made it here8" << std::endl;

      //if we don't have ks centers, then a center lost all its points, which is almost certainly an error
      //so we abort
      if (pts_grouped_by_center.size() != ks) {
        std::cout << "not the right size pgbc " << pts_grouped_by_center.size()  << std::endl;
        abort();
      }

      int debug_gbc = 0;
      for (size_t i = 0; i < pts_grouped_by_center.size(); i++) {
        debug_gbc += pts_grouped_by_center[i].second.size();
      }
      if (debug_gbc != ns) {
        std::cout << "pts grouped by center does not add to n, aborting " << std::endl;
  
      }

      //std::cout << "made it here9" << std::endl;

      //TODO only do update centers for unstable centers *mix with stable
      // UPDATING RADIUS

      parlay::parallel_for(0,ks,[&] (uint8_t i) {
        uint8_t picked_center_coord = pts_grouped_by_center[i].first;
        if (centers[picked_center_coord].stable == false) {
          auto list_of_dists = parlay::map(pts_grouped_by_center[i].second, [&] (int j) {
          return D.distance(v+j*ds,centers[picked_center_coord].tcoordinates.begin(),ds);
        });
        centers[picked_center_coord].radius = *parlay::max_element(list_of_dists);

        }
        
      } );

      //making all the radius super large, see if that fixes error (means we check all centers)
      //TODO remove this line
      // for (uint8_t i = 0; i < ks; i++) {
      //   centers[i].radius = std::numeric_limits<float>::max();
      // }

      
      // DEBUGGING
      //printing radius
      // std::cout << "printing radius" << std::endl;

      // for(int i = 0; i < ks; i++){
      //   std::cout << centers[i].radius << std::endl;
      // }

    //TODO make parallel later
    //TODO only compare for i < j 
    //TODO for iter  >= 1, use Alg 2 to reduce center-center distance calc
    //fill in the distance matrix
    for (uint8_t i = 0; i < ks; i++) {
      for (uint8_t j = 0; j < ks; j++) {
        dist_matrix[i][j] = D.distance(centers[i].coordinates.begin(),centers[j].coordinates.begin(),ds);
        //Debugging: print out the dist matrix
        //std::cout << "dist mat " << static_cast<int>(i) << " " << static_cast<int>(j) << ": " << dist_matrix[i][j] << std::endl;

      }
    }

  

    //std::cout << "made it here11" << std::endl;
    //set old neighbors to be the previous neighbors
    for (uint8_t i = 0; i < ks; i++) {
      //TODO is the below line a deep copy? or is a for loop needed
      //old_neighbors = neighbors; //copy here OK
      centers[i].old_neighbors = parlay::sequence<uint8_t>(centers[i].neighbors.size());
      for (uint8_t j = 0; j < centers[i].neighbors.size(); j++) {
        centers[i].old_neighbors[j] = centers[i].neighbors[j];
        
      }
      centers[i].neighbors = parlay::sequence<uint8_t>(0);
    }

    //std::cout << "made it here12" << std::endl;

    //init neighbor set
    for (uint8_t i = 0; i < ks; i++) {
      for (uint8_t j = 0; j < ks; j++) {
        if (i == j) {
          continue;
        }

        if (dist_matrix[i][j] < 2 * centers[i].radius) {
          centers[i].neighbors.push_back(j);
        }

      }
    }

    //std::cout << "made it here13" << std::endl;

    //sort neighbors 
    //TODO make sorting more efficient? (use parlay or a comparator?)
    for (uint8_t i = 0; i < ks; i++) {
     
      //DEBUGGING: Print out all of the neighbors of each center, and the radii
      // for (uint8_t j = 0; j < centers[i].neighbors.size(); j++) {
      //   std::cout << "(center, neighbor id, dist): " << "(" << static_cast<int>(i) << ", "<< static_cast<int>(centers[i].neighbors[j]) <<
      //   ", " << dist_matrix[i][centers[i].neighbors[j]] << ") " << std::endl;
      // }
      // std::cout << std::endl;

    
      sort_neighbor(centers,dist_matrix,i);
    }

    //debugging test
    //confirm neighbors are sorted
    /*
    for (uint8_t i = 0; i < ks; i++) {
      for (uint8_t j = 0; j < centers[i].neighbors.size(); j++) {
        std::cout << "(center, neighbor id, dist): " << "(" << static_cast<int>(i) << ", "<< static_cast<int>(centers[i].neighbors[j]) <<
        ", " << dist_matrix[i][centers[i].neighbors[j]] << ") " << std::endl;
      }
      std::cout << std::endl;

    }
    */
 
    //std::cout << "reached beginning of ball" << std::endl;
    //BALL ASSIGN

    //old assignments are the asg from the previous iter
    parlay::parallel_for(0, ns, [&](int i) {
        asgs_old[i] = asgs[i]; //set old the number from the previous iter
      
      });


    //TODO have array of tight bound from point to nearest center

    // to get where the center is stored
    //we want indices[i] to give the index in pts_grouped_by_center where
    //i is stored
      parlay::sequence<uint8_t> indices(ks);
      for(int i = 0; i < ks; i++){
        for(int j = 0; j < ks; j++){
            if(i == pts_grouped_by_center[j].first){
              indices[i] = j; 
              //std::cout << i << " corrosponds to " << j << std::endl;
              break;
            }
        }
      }
    //TODO look over the centers in parallel?
    //for each center
    //for (uint8_t i = 0; i < ks; i++) {
    parlay::parallel_for(0,ks,[&] (size_t i) {
      
      //help with debugging, make lines of code easier to read
   //   center& main_center = centers[i];
   //   parlay::sequence<uint8_t>& main_neighbors = centers[i].neighbors;

      bool flag_neighbors_stable = true;
      //TODO is iter >= 2 correct here?
      if (iter >= 2 && centers[i].neighbors == centers[i].old_neighbors && centers[i].stable) {
        for (uint8_t j = 0; j < centers[i].neighbors.size(); j++) {
          if (centers[centers[i].neighbors[j]].stable == false) {
            flag_neighbors_stable = false;
            break;
          }
        }
        
      }
      else {
        flag_neighbors_stable = false;
      }

      if (flag_neighbors_stable == false) {

        //find min dist to any other center
        //TODO use a map then min_element
       // float min_dist_to_other_neighbors = std::numeric_limits<float>::max();
       if (centers[i].neighbors.size() == 0) {
        std::cout << "no neighbors, flag should be stable, aborting" << std::endl;
        abort();
       }
        float min_dist_to_other_neighbors = dist_matrix[i][centers[i].neighbors[0]];
        // std::cout << "min dist to other neighbors: " <<  min_dist_to_other_neighbors << std::endl;
        //std::cout << "REACHED BALL ASSIGN" << std::endl;

        if (pts_grouped_by_center[indices[i]].first != i) {
          std::cout << "indices map wrong, aborting" << std::endl;
          abort();
        }

        //for each point belonging to that center
        parlay::sequence<int>& pts = pts_grouped_by_center[indices[i]].second;

        if (pts.size() != centers[i].old_num_members) {
          std::cout << "size old num members mismatch " << pts.size() << " " << centers[i].old_num_members << std::endl;
          abort();
        }
        //std::cout << "Center: " << static_cast<int>(i) << std::endl;
        /*
        for(int y = 0; y < pts.size(); y++){
            std::cout << pts[y] << std::endl;
        }
        */
      //  std::cout << pts_grouped_by_center[indices[i]].second.size() << ", " << centers[i].old_num_members << std::endl; 
      //for each point belonging to the center i, find its closest center
      parlay::parallel_for(0,pts.size(),[&] (size_t j) {
      //for (size_t j = 0; j < pts.size(); j++) {

            size_t ind = pts[j];
            
          //  std::cout << ind << std::endl;

            float dist_to_closest = D.distance(v+ind*ds,centers[i].tcoordinates.begin(),ds);
   
            if (dist_to_closest > 0.5 * min_dist_to_other_neighbors){
           
              float min_dist = dist_to_closest;
              uint8_t new_center = i;
              for(uint8_t u = 0; u < centers[i].neighbors.size(); u++){
           
                if (dist_to_closest < 0.5 * dist_matrix[i][centers[i].neighbors[u]]){
                  break;
                }
                if (min_dist > D.distance(v+ind*ds,centers[centers[i].neighbors[u]].tcoordinates.begin(),ds)){
                    min_dist = D.distance(v+ind*ds,centers[centers[i].neighbors[u]].tcoordinates.begin(),ds);
                    new_center = centers[i].neighbors[u];
                }
                
              }
              //std::cout << "min dist: " << min_dist << std::endl;
              asgs[ind] = new_center;
            }

        });
      
      }

    });

    //std::cout << "Finished BALL ASSIGN" << std::endl;



      /*
      //STEP 1: ASSIGN
      parlay::parallel_for(0, ns, [&](int i) {
        asgs_old[i] = asgs[i]; //set old the number from the previous iter
        auto distances = parlay::delayed::map(centers, [&] (center& q) {
          return D.distance(v+i*ds,q.tcoordinates.begin(),ds);
        });
        asgs[i] = min_element(distances) - distances.begin();
        
      });

      

      std::cout << "made it here 9 " << std::endl;
      */
     assign_time = t2.next_time();

      //STEP 2: UPDATE

      //STEP 2.1: PREPROCESSING

      //printing out all of asgs
      //specific for n=100
      // std::cout << "printing all of asgs" << std::endl;
      // for (int i = 0; i < 50; i++) {
      //   for (int j = 0; j < 20; j++) {
      //     std::cout << static_cast<int>(asgs[i*20 + j]) << " " ;

      //   }
      //   std::cout << std::endl;
      // }

      //use a histogram to get new num members info, this is faster
      auto new_centers_dist = histogram_by_key(parlay::map(rang,[&] (int i) {
        return asgs[i];
      }));

      if (new_centers_dist.size() != ks) {
      std::cout << "center not represented aborting " << new_centers_dist.size() << std::endl;
      abort();
    }


      for (uint8_t i = 0; i < ks; i++) {
        //LOL .size() to get the # of members!
        centers[new_centers_dist[i].first].new_num_members = new_centers_dist[i].second; //I think this is how it is used? 
        
        
      }

      //std::cout << "made it here1" << std::endl;

      //update whether or not each center has changed

      //TODO make more efficient? (if possible?) 
      parlay::parallel_for(0,ks,[&] (uint8_t i) {
        centers[i].has_changed = false;
      });

      //std::cout << "made it here 40" << std::endl;

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

      //std::cout << "made it here2" << std::endl;

      //STEP 2.2: Actually calculate new centers

      compute_centers_compare(v,ns,ds,ks,asgs,asgs_old,centers,center_calc_float,center_calc,rang);

      //STEP 2.3: Calculate info for next iter, copy over centers

      parlay::sequence<float> deltas = parlay::tabulate(ks,[&] (size_t i) {
        return D.distance(centers[i].coordinates.begin(),center_calc_float + i*ds,ds);
      });
      float total_diff = *parlay::max_element(deltas); //does max work like this?

      //copy back over centers
      parlay::parallel_for(0,ks,[&](size_t i) {
        for (unsigned short int j = 0; j < ds; j++) {
          centers[i].coordinates[j] = center_calc_float[i*ds+j];
          centers[i].tcoordinates[j] = centers[i].coordinates[j];
        }
        centers[i].delta = deltas[i];
      });


      for(int i = 0; i < ks; i++){
        
        centers[i].old_num_members = centers[i].new_num_members;
      }

      update_time=t2.next_time();


      //calculate msse
      //use tcoordinates! (so no casting happens inside)
      //FIXME when tcoordinates not used distances go to 0, 
      //bug in cross type code? 
      double msse = parlay::reduce(parlay::map(rang,[&] (int i) {
        return static_cast<double>(D.distance(v+i*ds,centers[asgs[i]].tcoordinates.begin(),ds));
      }))/ns;

      //std::cout << "msse is " << msse << std::endl;


      setup_time += t2.next_time();

      logger.add_iteration(assign_time,update_time,msse,0,0,deltas,setup_time);
      assign_time=0;
      update_time=0;
      setup_time=0;


      if (total_diff <= epsilon) {
        break;
      }

      // debug_num_members = 0;
      // debug_num_members_new = 0;
      // std::cout << "after iter "<< iter << std::endl;
      // for (uint8_t i = 0; i < ks; i++) {
      //   std::cout << "center " << static_cast<int>(i) << "has " << centers[i].old_num_members << std::endl;
      //   debug_num_members += centers[i].old_num_members;
      //   debug_num_members_new += centers[i].new_num_members;
      // }
      // if (debug_num_members != ns || debug_num_members_new != ns) {
      //   std::cout << "members don't add up, aborting" << std::endl;
      //   std::cout << debug_num_members << " " << debug_num_members_new << std::endl;
      //   abort();
      //  }



    }
    parlay::parallel_for(0,ns,[&] (int i) {
      asg[i] = asgs[i];
    });
    parlay::parallel_for(0,ks,[&] (size_t i) {
      for (unsigned short int coord = 0; coord < ds; coord++) {
        c[i*d+coord] = centers[i].coordinates[coord];
      }
    });


  }
};


#endif //SKLN