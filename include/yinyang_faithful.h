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

template<typename T>
struct YinyangFaithful {

    bool DEBUG_VD = false;

    struct point {

        parlay::slice<T*, T*> coordinates; // the coordinates of the point
        size_t best; // the index of the best center for the point
        size_t id;
        float ub;
        parlay::sequence<float> lb; //lower bounds from a point to each group
        float global_lb; //global lower bound

    
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

        float max_drift;
    };


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

    //TODO: If the type is float, then run a separate branch, 
    //no need to convert to float

    //find the distance from a point to each group (closest and 2nd closest)
    void distances_to_groups(const point&p, parlay::sequence<center>& centers,
    parlay::sequence<group> groups,
    Distance& D, size_t t, size_t group_size, 
    parlay::sequence<float>& distances, parlay::sequence<float> distances2nd) {

        //convert point coordinates to float buffer
        const int d = p.coordinates.size();
        float buf[d];
        T* it = p.coordinates.begin();
        for (int i = 0; i < d; i++) {
            buf[i]=*it;
            it += 1; //add 1 for next?
        }

        //for each group
        parallel_for(0,t,[&] (size_t i) {
            //calculate the distance from each group member to the point
            auto dist = parlay::map(groups[i].center_ids, 
            [&](size_t j) {
            return D.distance(buf, make_slice(centers[j].coordinates).begin(),d);
            });

            //find the closest and 2nd closest center
            //TODO make parallel
            size_t min_index = -1;
            size_t min2_index = -1; //2nd smallest
            if (dist[0] < dist[1]) {
                min_index = 0;
                min2_index = 1;
            }
            else {
                min_index = 1;
                min2_index = 0;
            }
            for (int j = 2; j < group_size; j++) {
                if (dist[j] < dist[min2_index]) {
                    if (dist[j] < dist[min_index]) {
                        min2_index = min_index;
                        min_index = j;
                    
                    }
                    else {
                        min2_index = j;
                    }
                 }

            }
            distances[i] = min_index;
            distances2nd[i] = min2_index;
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


        

        
        size_t group_size = 10;
        size_t t = k/group_size; //t is the number of groups
        if (group_size * t != k) {
            std::cout << "not exactly divisible t k, aborting" << std::endl;
            abort();
        }
        if (group_size <= 1) {
            std::cout << "Group size 1 not supported right now, aborting\n";
            abort();
        }


        //cluster on the groups initially using NaiveKmeans
        float* group_centers = new float[t * d];
        size_t* group_asg = new size_t[k];
        LazyStart<float> init;
        init.initialize(c,k,d,t,group_centers,group_asg,D);
        NaiveKmeans<float> run;
        run.cluster(c,k,d,t,
            group_centers, group_asg,D, 5, 0.0001);

        parallel_for(0,k,[&] (size_t i) {
            centers[i].group_id = group_asg[i];
        });

        parlay::sequence<group> groups(t);


        //Sadly sequential assigning the groups
        for (int i = 0; i < k; i++) {
            groups[centers[i].groups_id].center_ids.push_back(i);
        }

        //confirm group assignment happened properly
        //(checking a necessary not sufficient condition)
        for (int i =0 ;i < t; i++) {
            if (groups[i].center_ids.size() != group_size) {
                std::cout << 
                "Group assignment went wrong, group is wrong size"
                << std::endl;
                abort();
            }
        }


        //Yinyang first does a normal iteration of kmeans:   
        
        //Assignment
        //distance from each group to a point
        parlay::sequence<parlay::sequence<float>> distances(n,parlay::sequence<float>(t));
        //distance from 2nd closest point in each group to a point
        parlay::sequence<parlay::sequence<float>> distances2nd(n,parlay::sequence<float>(t));
        // Assign each point to the closest center
        parlay::parallel_for(0, pts.size(), [&](size_t i) {
            distances_to_groups(pts[i],centers,groups,D,t,group_size,distances,
            distances2nd);

            pts[i].lb = parlay::sequence<float>(t); //initialize the lower bound sequence
        });

        //actually take best
        parlay::parallel_for(0,n,[&] (size_t i) {
            pts[i].best=min_element(distances[i])-distances[i].begin();
            pts[i].ub = distances[pts[i].best];
            //could do in parallel? TODO
            for (int j = 0; j < t; j++) {
                //in the group containing the best, you take the 2nd best from
                //the group
                if (pts[i].best==j) {
                    pts[i].lb[j] = distances2nd[i][j];
                }
                else {
                    pts[i].lb[j] = distances[i][j];
                }
                
            }
        });

       

        size_t iters = 1;
        
        //Step 3: Repeat until convergence
        while (true) {

             // Compute new centers
            parlay::sequence<center> new_centers = compute_centers_vd(pts, n, d, k, 
                centers);

            // Check convergence
            total_diff = 0.0;
            for (size_t i = 0; i < k; i++) { 
            new_centers[i].delta = D.distance(centers[i].coordinates.begin(), 
            new_centers[i].coordinates.begin(),k);
            total_diff += new_centers[i].delta;
            }

            centers = std::move(new_centers);

            //for each group, get max drift for group
            parlay::parallel_for(0,t,[&] (size_t i) {
                auto drifts = map(groups[i].center_ids, [&] (size_t j) {
                    return centers[j].delta;

                });
                groups[i].max_drift = max_element(drifts);

            });


            //3.2: Group filtering

            //update bounds
            parlay::parallel_for(0,n,[&](size_t i) {
                pts[i].ub += centers[pts[i].best].delta; 
                pts[i].global_lb = max(0,pts[i].lb-groups[j].max_drift);
                for (int j = 0; j < t; j++) {
                    pts[i].lb[j] = max(0,pts[i].lb-groups[j].max_drift);
                    if (pts[i].global_lb > pts[i].lb[j]) {
                        pts[i].global_lb = pts[i].lb[j];
                    } //continue here TODO
                }


            });








            //convergence check
            if (iters >= max_iter || total_diff < epsilon) break;
            iters += 1;


        }
        

        delete[] group_centers; //memory cleanup
        delete[] group_asg;

    }
    
};