#include <tuple>
#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/delayed.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"
#include "NSGDist.h"

template<typename T>
struct YinyangFaithful {

    struct point {

        parlay::slice<T*, T*> coordinates; // the coordinates of the point
        size_t best; // the index of the best center for the point
        size_t id;
        float ub;
        float lb;

    
        point() : best(-1), coordinates(nullptr, nullptr), 
        ub(std::numeric_limits<float>::max()), lb(-1) {
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates) : best(chosen), 
        coordinates(coordinates.begin(),coordinates.end()), 
        ub(std::numeric_limits<float>::max()), lb(-1) {
        
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates, float ub, 
        float lb) : best(chosen), 
        coordinates(coordinates.begin(),coordinates.end()), 
        ub(ub), lb(lb) {
        
        }



    };

    struct center {
        size_t id; // a unique (hopefully) identifier for the center
        size_t group_id; //the id of the group that the center belongs to
        parlay::sequence<float> coordinates; // the pointer to coordinates of the center
        float delta;

        center(size_t id, parlay::sequence<float> coordinates) : id(id) {
        
            this->coordinates = coordinates;
        }

        center() : id(-1) {
        
        }

    };

    struct group {
        size_t id;
        //store the ids of all the centers 
        //belonging to this group
        parlay::sequence<size_t> center_ids;
    }


    //We copy in the naive kmeans code because Yinyang first does a normal 
    //iteration of kmeans

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
    for (int i = 0; i < k; i++) {
        new_centers[i].id = i;
        new_centers[i].coordinates=parlay::sequence<T>(d,4);
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

    //distances to each center
    //distances put in the distances argument, a sequence of size k
    void distances_to_centers(const point& p, parlay::sequence<center>& 
    centers, Distance& D, parlay::sequence<float>& distances) {
        size_t d = p.coordinates.size();
        parallel_for(0,k,[&] (size_t i) {
            distances[i] = D.distance(p.coordinates.begin(),
            make_slice(q.coordinates).begin(),d);
        });

    }



    void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
    Distance& D, size_t max_iter, double epsilon) {


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

        size_t t = k/10; //t is the number of groups

        //cluster on the groups initially using NaiveKmeans
        float* group_centers = new float[t * d];
        size_t* group_asg = new size_t[k];

        Kmeans<float,LazyStart<float>,NaiveKmeans<float>>(c,k,d,t,
            group_centers, *D, 5, 0.0001);

        parallel_for(0,k,(size_t i) {
            centers[i].group_id = group_asg[i];
        })

        //Yinyang first does a normal iteration of kmeans:   
        
        //Assignment
        parlay::sequence<float> distances(k);
        distances_to_centers()
        // Assign each point to the closest center
        parlay::parallel_for(0, pts.size(), [&](size_t i) {
        pts[i].best = closest_point_vd(pts[i], centers,D);
        });

        // Compute new centers
        parlay::sequence<center> new_centers = compute_centers_vd(pts, n, d, k, 
            centers);

        // Check convergence
        total_diff = 0.0;
        for (size_t i = 0; i < k; i++) { 
        float diff = D.distance(centers[i].coordinates.begin(), 
        new_centers[i].coordinates.begin(),k);
        total_diff += diff;
        }

        centers = std::move(new_centers);

        delete[] group_centers; //memory cleanup
        delete[] group_asg;

    }
    
};