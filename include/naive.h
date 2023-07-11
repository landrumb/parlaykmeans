/* 
    naive implementation of Lloyd's algorithm for k-means clustering
 */

#ifndef NAIVE
#define NAIVE

#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/delayed.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"
#include "utils/NSGDist.h"
#include "utils/kmeans_bench.h"

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
    bool DEBUG_VD = false;

    struct point {

        size_t best; // the index of the best center for the point
        parlay::slice<T*, T*> coordinates; // the coordinates of the point

    
        point() : best(-1), coordinates(nullptr, nullptr) {
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates) : best(chosen), 
        coordinates(coordinates.begin(),coordinates.end()) {
        
        }
    };

    struct center {
        size_t id; // a unique (hopefully) identifier for the center
        parlay::sequence<float> coordinates; // the pointer to coordinates of the center
    
        center(size_t id, parlay::sequence<float> coordinates) : id(id) {
        
            this->coordinates = coordinates;
        }

        center() : id(-1) {
        
        }

    };


    //useful helper functions
    //corresponding prints

    void print_point(const point& p) {
      if constexpr(std::is_same<T,uint8_t>()==true) {
        std::cout << "Po: " << p.best << std::endl;
        for (T* i = p.coordinates.begin(); i != p.coordinates.end(); i = std::next(i) ) {
          std::cout << static_cast<int>(*i) << " ";
        }
        std::cout << std::endl;

      }
      else {
        std::cout << "Po: " << p.best << std::endl;
        for (T* i = p.coordinates.begin(); i != p.coordinates.end(); i = std::next(i) ) {
          std::cout << *i << " ";
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
      
    }


    //MUST PASS DISTANCE BY REFERENCE NOT COPY
    //put the coordinates of p onto the stack (in buf) for the calculation
    size_t closest_point_vd(const point& p, 
    parlay::sequence<center>& centers, Distance& D) {
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
       
        const int d = p.coordinates.size();
        float buf[2048];
        T* it = p.coordinates.begin();
        for (int i = 0; i < d; i++) {
            buf[i]=*it;
            it += 1; //add 1 for next?
        }
     
        auto distances = parlay::delayed::map(centers, [&](center& q) {
            return D.distance(buf, make_slice(q.coordinates).begin(),d);
        });

        //C++ won't auto-cast
        // auto distances = parlay::delayed::map(centers, [&](center& q) {
        //     return D.distance(p.coordinates.begin(), make_slice(q.coordinates).begin(),d);
        // });

        if (DEBUG_VD) {
          std::cout << "distance printing" << std::endl;
          for (size_t i = 0; i < distances.size(); i++) {
            std::cout << distances[i] << " " ;
          }
          std::cout << std::endl;
        }

        return min_element(distances) - distances.begin();

    }

  }


//compute centers calculates the new centers
//returns: a sequence of centers
parlay::sequence<center> compute_centers_ec(const 
parlay::sequence<point>& pts, size_t n, size_t d, size_t k, const parlay::sequence<center>& centers) {

    parlay::sequence<parlay::sequence<size_t>> indices(k);
    
    parlay::sequence<center> new_centers(k);
    for (size_t i = 0; i < k; i++) {
        new_centers[i].id = i;
        new_centers[i].coordinates=parlay::sequence<T>(d,4);
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
parlay::sequence<center> compute_centers_vd(
  const parlay::sequence<point>& pts, size_t n, size_t d, size_t k, 
  const parlay::sequence<center>& centers) {

    if constexpr(std::is_same<T,float>() == true) {
        //run the normal compute centers
        return compute_centers_ec(pts,n,d,k,centers);
    }

    parlay::sequence<parlay::sequence<size_t>> indices(k);
    
    parlay::sequence<center> new_centers(k);
    for (size_t i = 0; i < k; i++) {
        new_centers[i].id = i;
        new_centers[i].coordinates=parlay::sequence<float>(d,4);
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

void cluster(T* v, size_t n, size_t d, size_t k, 
float* c, size_t* asg, Distance& D,  size_t max_iter, double epsilon) {

  //format the data according to our naive run
  parlay::sequence<point> pts = parlay::tabulate<point>(n, [&] (size_t i) {
      return point(asg[i],parlay::slice(v+i*d, v+i*d + d));

  });

  //create the centers
  parlay::sequence<center> centers = parlay::tabulate<center>(k, [&] (size_t i) {
    return center(i, parlay::sequence<float>(d));
  });

  //fill in the centers
  parlay::parallel_for(0,k, [&] (size_t i) {
    for (size_t j = 0; j < d; j++) {
      centers[i].coordinates[j] = *(c + i*d + j);
    }

  });

  //the actual naive run
  kmeans_vd(pts,n,d,k,centers,D,max_iter,epsilon);

  //put our data back 
  parlay::parallel_for(0,k,[&] (size_t i) {
    for (size_t j = 0; j < d; j++) {
      c[i*d + j] = centers[i].coordinates[j];
    }
  });
  parlay::parallel_for(0,n,[&] (size_t i) {
      asg[i] = pts[i].best;
  });

}


double kmeans_vd(parlay::sequence<point>& pts, size_t n, size_t d, size_t k, 
parlay::sequence<center>& centers, Distance& D, size_t max_iterations, 
double epsilon)
{

    std::cout << "running vd" << std::endl;


    if (d > 2048) {
      std::cout << "d greater than 2048, too big, printing d: " << 
      d << std::endl;
      abort();
    }

    parlay::internal::timer t = parlay::internal::timer();
    t.start();
    
     if(pts.size() == 0){ //default case
        return -1;
    }


  
    parlay::sequence<float> buf1(50,1);
    parlay::sequence<float> buf2(50,2);
    std::cout << "dist " << D.distance(make_slice(buf1).begin(),
    make_slice(buf2).begin(),50) << std::endl;
    std::cout << "finished with dist\n";
  size_t iterations = 0;

  float total_diff = 0;

  // t.next("finished init");

  while (iterations < max_iterations) {

    if (DEBUG_VD) {
         std::cout << "centers: " << iterations << std::endl;
    for (size_t i = 0; i < k; i++) {
       print_center(centers[i]);        
    }

    }
    

   
    // std::cout << "iter" << iterations << std::endl;
    iterations++;

    if (DEBUG_VD) {
        std::cout << "center printing" << std::endl;
        for (size_t i = 0; i < k; i++) {
            print_center(centers[i]);
        }
        std::cout << std::endl << std::endl;

    }
     
    t.next("About to do closest points");
    // Assign each point to the closest center
    parlay::parallel_for(0, pts.size(), [&](size_t i) {
      pts[i].best = closest_point_vd(pts[i], centers,D);
    });

    t.next("Finished closest points");

    //parallel version for debugging:
    // for (size_t i = 0; i < pts.size(); i++) {
    //   pts[i].best = closest_point_vd_seq(pts[i], centers,*D);

    // }

    std::cout << "past closest pt\n";

    // Compute new centers
    parlay::sequence<center> new_centers = compute_centers_vd(pts, n, d, k, centers);
    t.next("Computed the centers");
    // Check convergence
    total_diff = 0.0;
    for (size_t i = 0; i < k; i++) { 
      float diff = D.distance(centers[i].coordinates.begin(), new_centers[i].coordinates.begin(),k);
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
  for (size_t i = 0; i < k; i++) {
    print_center(centers[i]);
  }
  std::cout << "Error" << total_diff << std::endl;
  std::cout << std::endl << std::endl;

    }
  

  return t.total_time(); //std::make_pair(centers,timer.total_time());

    //return std::make_pair(centers,timer.total_time());

}


};

template <class T>
struct Naive {
  void operator()(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg,
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
    parlay::internal::timer t = parlay::internal::timer();
    t.start();

    float* temp_centers = new float[k * d];
    threadlocal::accumulator<double> squared_errors = threadlocal::accumulator<double>();

    for (size_t i = 0; i < max_iter; i++) {
      // Assign each point to the closest center
      parlay::parallel_for(0, n, [&](size_t i) {
        // float min_dist = std::numeric_limits<float>::max();
        float min_dist = D.distance(v + i * d, c + asg[i] * d, d);
        for (size_t j = 0; j < k; j++) {
          float dist = D.distance(v + i * d, c + j * d, d);
          if (dist < min_dist) {
            min_dist = dist;
            asg[i] = j;
          }
        }
        squared_errors.add(min_dist);
      });
      
      float assignment_time = t.next_time();

      // Update centers
      threadlocal::accumulator<double>* new_centers = new threadlocal::accumulator<double>[k * d];
      threadlocal::accumulator<size_t>* assignments = new threadlocal::accumulator<size_t>[k];
      parlay::parallel_for(0, n, [&](size_t i) {
        for (size_t j = 0; j < d; j++) {
          new_centers[asg[i] * d + j].add(v[i * d + j]);
          assignments[asg[i]].increment();
        }
      });
      for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < d; j++) {
          temp_centers[i * d + j] = static_cast<float>(new_centers[i * d + j].total() / assignments[i].total());
        }
      }
      
      auto deltas = parlay::tabulate(k, [&](size_t i) {
        return D.distance(temp_centers + i * d, c + i * d, d);
      });

      // copy temp centers into c
      std::copy(temp_centers, temp_centers + k * d, c);

      float max_delta = parlay::reduce(deltas, parlay::maxm<float>());
      float update_time = t.next_time();
      logger.add_iteration(assignment_time, update_time, squared_errors.total() / n, 0, 0, deltas);
      if (max_delta <= epsilon) {
        break;
      }

      // reset accumulators
      for (size_t i = 0; i < k * d; i++) {
        new_centers[i].reset();
      }
      for (size_t i = 0; i < k; i++) {
        assignments[i].reset();
      }
      squared_errors.reset();
    }

    // delete[] temp_centers;
    // delete[] new_centers;
    // delete[] assignments;

    return;
  };

  std::string name() { return "naive"; };

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg,
  Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
    this->operator()(v, n, d, k, c, asg, D, logger, max_iter, epsilon);
  };
};


#endif //NAIVE