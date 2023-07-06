/* 
    naive implementation of Lloyd's algorithm for k-means clustering
 */

#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/delayed.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

template <typename T>
struct NaiveKmeans {

    //print out a ton of debug info if true
    private: bool DEBUG_VD = false;

    private: struct point {

        parlay::slice<T*, T*> coordinates; // the coordinates of the point
        size_t best; // the index of the best center for the point

    
        point() : best(-1), coordinates(nullptr, nullptr) {
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates) : best(chosen) {
          this->coordinates=coordinates;

        }
    };

    private: struct center {
        size_t id; // a unique (hopefully) identifier for the center
        parlay::sequence<float> coordinates; // the pointer to coordinates of the center
    
        center(size_t id, parlay::sequence<float> coordinates) : id(id) {
        
            this->coordinates = coordinates;
        }

        center() : id(-1) {
        
        }

    };

    //MUST PASS DISTANCE BY REFERENCE NOT COPY
    //put the coordinates of p onto the stack (in buf) for the calculation
    private: size_t closest_point_vd(const point<T>& p, sequence<center>& 
    centers, Distance& D) {
    if constexpr(std::is_same<T,float>() == true) {
     

          int d = p.coordinates.size();

        //no need to convert with a buffer
        auto distances = parlay::delayed::map(centers, [&](center& q) {
            return D.distance(p.coordinates.begin(), 
              make_slice(q.coordinates).begin(),d);
        });

        return min_element(distances) - distances.begin();

        
    }
    else {
       
        int d = p.coordinates.size();
        float buf[d];
        T* it = p.coordinates.begin();
        for (int i = 0; i < d; i++) {
            buf[i]=*it;
            it += 1; //add 1 for next?
        }
     
        auto distances = parlay::delayed::map(centers, [&](center& q) {
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


//compute centers calculates the new centers
//returns: a sequence of centers
template<typename T> sequence<center<T>> compute_centers_ec(const 
sequence<point<T>>& pts, size_t n, size_t d, size_t k, const sequence<center<T>>& centers) {

    sequence<sequence<size_t>> indices(k);
    
    sequence<center<T>> new_centers(k);
    for (int i = 0; i < k; i++) {
        new_centers[i].id = i;
        new_centers[i].coordinates=sequence<T>(d,4);
    }

    for (int i = 0; i < n; i++) {
            
        indices[pts[i].best].push_back(i); //it's called best not id!!!

    }
    if (DEBUG_ED) {
        std::cout << "Debugging: printing out center counts:\n";
        for (int i = 0; i < k; i++) {
        std::cout << indices[i].size() << std::endl;
        }

    }
    
    parlay::parallel_for (0, k*d, [&] (size_t icoord){
        size_t i = icoord / d;
        size_t coord = icoord % d;

        //std::cout<<"icoord " << icoord << "i : " << i << "coord: " << 
        //  coord << std::endl;
        //new_centers[i].coordinates[coord] = 
        //  anti_overflow_avg(map(new_centers[i].points,[&] 
        //  (size_t ind) {return pts[ind].coordinates[coord];}  ));

        //if there are no values in a certain center, just keep the center
        // where it is
        if (indices[i].size() > 0) { //anti_overflow_avg or reduce??
            new_centers[i].coordinates[coord] = reduce(map(indices[i],[&] 
            (size_t ind) {
              return pts[ind].coordinates[coord];})) 
              / indices[i].size(); //normal averaging now

        }
        else { 
            new_centers[i].coordinates[coord] = centers[i].coordinates[coord];
        }

    });


    return new_centers;


}

//compute centers calculates the new centers
//returns: a sequence of centers
template<typename T> sequence<center> compute_centers_vd(
  const sequence<point<T>>& pts, size_t n, size_t d, size_t k, const sequence<center>& centers, 
  ) {

    if constexpr(std::is_same<T,float>() == true) {
        //run the normal compute centers
        return compute_centers_ec(pts,n,d,k,centers);
    }

    sequence<sequence<size_t>> indices(k);
    
    sequence<center> new_centers(k);
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

      
        //if there are no values in a certain center, just keep the center 
        //where it is
        if (indices[i].size() > 0) { //anti_overflow_avg or reduce?? 
        //note the static cast to double here, because points are whatever
            new_centers[i].coordinates[coord] = reduce(map(indices[i],[&] 
            (size_t ind) {return static_cast<float>(
              pts[ind].coordinates[coord]);})) / indices[i].size(); 
              //normal averaging now

        }
        else { 
            new_centers[i].coordinates[coord] = centers[i].coordinates[coord];
        }


    });


    return new_centers;


}

template <typename T> void cluster(T* v, size_t n, size_t d, size_t k, 
float* c, size_t* asg, Distance& D, size_t max_iter, double epsilon) {

  //format the data according to our naive run
  parlay::sequence<point<T>> pts = parlay::tabulate<point<T>>(n, [&] (size_t i) {
      return point<T>(asg[i],parlay::slice(v+i*d, v+i*d + d));

  })

  //create the centers
  parlay::sequence<center> centers = parlay::tabulate<center>(k, [&] (size_t i) {
    return center(i, parlay::sequence<float>(d));
  });

  //fill in the centers
  parlay::parallel_for(0,k, [&] (size_t i) {
    for (int j = 0; j < d; j++) {
      centers[i].coordinates[j] = *(c + i*d + j);
    }

  });

  //the actual naive run
  kmeans_vd(pts,n,d,k,centers,D,max_iterations,epsilon);

  //put our data back 
  parlay::parallel_for(0,k,[&] (size_t i)) {
    for (int j = 0; j < d; j++) {
      c[i*d + j] = centers[i].coordinates[j];
    }
  }
  parlay::parallel_for(0,n,[&] (size_t i) {
      asg[i] = pts[i].best;
  });

}




double kmeans_vd(parlay::sequence<point<T>>& pts, size_t n, size_t d, size_t k, 
parlay::sequence<center>& centers, Distance& D, size_t max_iterations, double epsilon)
{

    std::cout << "running vd" << std::endl;

    parlay::internal::timer t = parlay::internal::timer();
    t.start();
    
     if(pts.size() == 0){ //default case
        return -1;
    }


  
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
    sequence<center> new_centers = compute_centers_vd(pts, n, d, k, centers);
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


}

