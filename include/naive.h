/* 
    naive implementation of Lloyd's algorithm for k-means clustering
 */


//print out a ton of debug info if true
bool DEBUG_VD = false;

#include "../parlay/random.h"
#include "../parlay/parallel.h"
#include "../parlay/primitives.h"
#include "../parlay/sequence.h"
#include "../parlay/slice.h"
#include "../parlay/delayed.h"
#include "../parlay/io.h"
#include "../parlay/internal/get_time.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <atomic>
#include <mutex>

template <typename T>
struct point {

    parlay::slice<T*, T*> coordinates; // the coordinates of the point
    size_t best; // the index of the best center for the point

    //comment out if needed
    float lb;
    float ub;

   
    point() : best(-1), coordinates(nullptr, nullptr) {
    }


};

template <typename T>
struct center {
    size_t id; // a unique (hopefully) identifier for the center
    parlay::sequence<T> coordinates; // the pointer to coordinates of the center
   
    center(size_t id, parlay::sequence<T> coordinates) : id(id) {
      
        this->coordinates = coordinates;
    }

    center() : id(-1) {
       
    }

    float change;

};

#include "types.h"
#include "center_creation.h"
#include "NSGDist.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <atomic>
#include <mutex>

using namespace parlay;


//MUST PASS DISTANCE BY REFERENCE NOT COPY
//put the coordinates of p onto the stack (in buf) for the calculation
template<typename T> size_t closest_point_vd(const point<T>& p, sequence<center<float>>& centers, Distance& D) {
    if constexpr(std::is_same<T,float>() == true) {
      assert(centers.size() > 0);
        assert(p.coordinates.size() == centers[0].coordinates.size());

          int d = p.coordinates.size();

        //no need to convert with a buffer
        auto distances = parlay::delayed::map(centers, [&](center<float>& q) {
            return D.distance(p.coordinates.begin(), make_slice(q.coordinates).begin(),d);
        });

        return min_element(distances) - distances.begin();

        
    }
    else {
        assert(centers.size() > 0);
        assert(p.coordinates.size() == centers[0].coordinates.size());
        int d = p.coordinates.size();
        float buf[d];
        T* it = p.coordinates.begin();
        for (int i = 0; i < d; i++) {
            buf[i]=*it;
            it += 1; //add 1 for next?
        }
     
        auto distances = parlay::delayed::map(centers, [&](center<float>& q) {
            return D.distance(buf, make_slice(q.coordinates).begin(),d);
        });

        if (DEBUG_VD) {
          std::cout << "distance printing" << std::endl;
          for (int i = 0; i < distances.size(); i++) {
            std::cout << distances[i] << " " ;
          }
          std::cout << std::endl;
        }
        return min_element(distances) - distances.begin();

    }

}


//MUST PASS DISTANCE BY REFERENCE NOT COPY
//put the coordinates of p onto the stack (in buf) for the calculation
//this version will try to be more parallel by doing a coordinate on every center then adding
//nope not faster (at least on this computer my Macbook)
template<typename T> size_t closest_point_vd_spanned(const point<T>& p, sequence<center<float>>& centers, Distance& D) {
    if constexpr(std::is_same<T,float>() == true) { //not changing this piece
      assert(centers.size() > 0);
        assert(p.coordinates.size() == centers[0].coordinates.size());

          int d = p.coordinates.size();

        //no need to convert with a buffer
        auto distances = parlay::delayed::map(centers, [&](center<float>& q) {
            return D.distance(p.coordinates.begin(), make_slice(q.coordinates).begin(),d);
        });

        return min_element(distances) - distances.begin();

        
    }
    else {
        int k = centers.size();
        int d = p.coordinates.size();

        sequence<sequence<float>> distances_expanded(k,sequence<float>(d));
        parallel_for(0,k*d,[&] (size_t icoord){
          int coord = icoord % d;
          int cen = icoord/d;
          float dist = static_cast<float>(p.coordinates[coord])-centers[cen].coordinates[coord];
          distances_expanded[cen][coord] = dist*dist;
        });

        auto  distances = parlay::delayed::map(distances_expanded, [&] (sequence<float> seq) {
          return reduce(seq);

        });

        if (DEBUG_VD) {
          std::cout << "distance printing" << std::endl;
          for (int i = 0; i < distances.size(); i++) {
            std::cout << distances[i] << " " ;
          }
          std::cout << std::endl;
        }
        return min_element(distances) - distances.begin();

    }

}





//put the coordinates of p onto the stack (in buf) for the calculation
//purely for debugging purposes, a sequential implementation easier to debug
template<typename T> size_t closest_point_vd_seq(const point<T>& p, sequence<center<float>>& centers, Distance& D) {
    if constexpr(std::is_same<T,float>() == true) {
      assert(centers.size() > 0);
        assert(p.coordinates.size() == centers[0].coordinates.size());

          int d = p.coordinates.size();

        //no need to convert with a buffer
        auto distances = parlay::delayed::map(centers, [&](center<float>& q) {
            return D.distance(p.coordinates.begin(), make_slice(q.coordinates).begin(),d);
        });

        return min_element(distances) - distances.begin();

        
    }
    else {
        assert(centers.size() > 0);
        assert(p.coordinates.size() == centers[0].coordinates.size());

        sequence<float> buf1(50,1);
        sequence<float> buf2(50,2);
        std::cout << "dist " << D.distance(make_slice(buf1).begin(),make_slice(buf2).begin(),50) << std::endl;
        std::cout << "finished with dist" << std::endl;
        int d = p.coordinates.size();
        float buf[d];
        T* it = p.coordinates.begin();
        for (int i = 0; i < d; i++) {
            buf[i]=*it;
            it += 1; //add 1 for next?
        }
        
        sequence<float> distances(centers.size());
        sequence<float> distances2(centers.size());
        for (int i = 0; i < centers.size(); i++) {
          // std::cout << "printing out first two points of buffer " << i << std::endl;
          // std::cout << buf[0] << " " << buf[1] << std::endl;
          std::cout << "about to call dist" << std::endl;
          distances[i] = D.distance(buf,make_slice(centers[i].coordinates).begin(),d);
          std::cout << "just called dist" << std::endl;
          distances2[i] = euclidean_squared(buf,centers[i].coordinates);
        }
        std::cout << "printing different distance calculations: " << std::endl;
        print_seq(distances);
        print_seq(distances2);
     
        assert(distances==distances2);
        assert(distances!=distances2);
        assert(distances[0]==distances2[0]);
     
        // auto distances = parlay::delayed::map(centers, [&](center<float>& q) {
        //     return D.distance(buf, make_slice(q.coordinates).begin(),d);
        // });

        if (DEBUG_VD) {
          std::cout << "distance printing" << std::endl;
          for (int i = 0; i < distances.size(); i++) {
            std::cout << distances[i] << " " ;
          }
          std::cout << std::endl;
        }
        return min_element(distances) - distances.begin();

    }

}


//compute centers calculates the new centers
//returns: a sequence of centers
template<typename T> sequence<center<float>> compute_centers_vd(const sequence<point<T>>& pts, const sequence<center<float>>& centers, size_t k, size_t d, size_t n) {

    if constexpr(std::is_same<T,float>() == true) {
        //run the normal compute centers
        return compute_centers_ec(pts,centers,k,d,n);
    }

    sequence<sequence<size_t>> indices(k);
    
    sequence<center<float>> new_centers(k);
    for (int i = 0; i < k; i++) {
        new_centers[i].id = i;
        new_centers[i].coordinates=sequence<float>(d,4);
    }


    for (int i = 0; i < n; i++) {
            
        indices[pts[i].best].push_back(i); //it's called best not id!!!

    }
    if (DEBUG_VD) {
        std::cout << "Debugging: printing out center counts:\n";
        for (int i = 0; i < k; i++) {
        std::cout << indices[i].size() << std::endl;
        }

    }
    
    
    parlay::parallel_for (0, k*d, [&] (size_t icoord){
        size_t i = icoord / d;
        size_t coord = icoord % d;

      
        //if there are no values in a certain center, just keep the center where it is
        if (indices[i].size() > 0) { //anti_overflow_avg or reduce?? //note the static cast to double here, because points are whatever
            new_centers[i].coordinates[coord] = reduce(map(indices[i],[&] (size_t ind) {return static_cast<float>(pts[ind].coordinates[coord]);})) / indices[i].size(); //normal averaging now

        }
        else { 
            new_centers[i].coordinates[coord] = centers[i].coordinates[coord];
        }


    });


    return new_centers;


}





template <typename T> double kmeans_vd(parlay::sequence<point<T>>& pts, parlay::sequence<center<float>>& centers, size_t k, size_t max_iterations, double epsilon){

    std::cout << "running vd" << std::endl;

    parlay::internal::timer t = parlay::internal::timer();
    t.start();
    

     if(pts.size() == 0){ //default case
        return -1;
    }

    size_t d = pts[0].coordinates.size();
    size_t n = pts.size();

    assert (k >= 32);
    assert (d >= 32);

    // if (d <= 32 || k <= 32) { //if d is small, vectorized code doesn't work so run a normal kmeans implementation
    //     std::cout << "using another type of kmeans run, not vectorization, as the data types aren't applicable" << std::endl;
    //     //copy the centers to doubles
    //     //need to properly copy the centers because we are passing by ref
    //     //FIXME because of the pass by reference the centers (floats) will look unchanged, even though the return value will be correct**
    //   parlay::sequence<center<double>> centers2(k,center<double>()); 
    //   for (int i = 0; i < k; i++) {
    //       centers2[i].id = centers[i].id;
    //       centers2[i].coordinates=sequence<double>(d);
    //       for (int j = 0; j < d; j++) {
    //           centers2[i].coordinates[j] = static_cast<double>(centers[i].coordinates[j]);
    //       }
    //   }

    //     return kmeans_double_center(pts,centers2,k,d,n).second;

    // }
    
    //if our data is in a double format, just run normal code, don't need a vectorized version
    if constexpr(std::is_same<T,double>() == true) {

       parlay::sequence<center<double>> centers2(k,center<double>()); 
        for (int i = 0; i < k; i++) {
            centers2[i].id = centers[i].id;
            centers2[i].coordinates=sequence<double>(d);
            for (int j = 0; j < d; j++) {
                centers2[i].coordinates[j] = static_cast<double>(centers[i].coordinates[j]);
            }
        }
        return kmeans_euc_dist(pts,centers2,k,d,n);
    }

    Distance* D = new Euclidian_Distance();

   
        sequence<float> buf1(50,1);
        sequence<float> buf2(50,2);
        std::cout << "dist " << (*D).distance(make_slice(buf1).begin(),make_slice(buf2).begin(),50) << std::endl;
        std::cout << "finished with dist\n";
  int iterations = 0;

  float total_diff = 0;

  t.next("finished init");

  while (iterations < max_iterations) {

    if (DEBUG_VD) {
         std::cout << "centers: " << iterations << std::endl;
    for (int i = 0; i < k; i++) {
       print_center(centers[i]);        
    }

    }
    

   
    std::cout << "iter" << iterations << std::endl;
    iterations++;

    if (DEBUG_VD) {
        std::cout << "center printing" << std::endl;
        for (int i = 0; i < k; i++) {
            print_center(centers[i]);
        }
        std::cout << std::endl << std::endl;

    }
     
    t.next("About to do closest points");
    // Assign each point to the closest center
    parlay::parallel_for(0, pts.size(), [&](size_t i) {
      pts[i].best = closest_point_vd(pts[i], centers,*D);
    });

    t.next("Finished closest points");

    //parallel version for debugging:
    // for (size_t i = 0; i < pts.size(); i++) {
    //   pts[i].best = closest_point_vd_seq(pts[i], centers,*D);

    // }

    std::cout << "past closest pt\n";

    // Compute new centers
    sequence<center<float>> new_centers = compute_centers_vd(pts, centers, k, d, n);
    t.next("Computed the centers");
    // Check convergence
    total_diff = 0.0;
    for (int i = 0; i < k; i++) { 
      float diff = D->distance(centers[i].coordinates.begin(), new_centers[i].coordinates.begin(),k);
      total_diff += diff;
    }

    centers = std::move(new_centers);

    std::cout << "difs " << total_diff << " " << epsilon << std::endl;

    if (total_diff <= epsilon) {
      break;
    }
    
  }
    if (DEBUG_VD) {
         std::cout << "center printing" << std::endl;
  for (int i = 0; i < k; i++) {
    print_center(centers[i]);
  }
  std::cout << "Error" << total_diff << std::endl;
  std::cout << std::endl << std::endl;

    }
  

  return t.total_time(); //std::make_pair(centers,timer.total_time());

    //return std::make_pair(centers,timer.total_time());

}