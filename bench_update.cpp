//file to rate the quality of different update methods
//goal to show that updates at high n are a problem
//for k-means and need to be addressed

#ifndef UPDATE_H
#define UPDATE_H

#include "include/utils/parse_command_line.h"
#include "include/utils/parse_files.h"
#include "include/utils/NSGDist.h"
#include "include/utils/kmeans_bench.h"

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/random.h"
#include "parlay/internal/get_time.h"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <set>
#include <atomic>
#include <utility>
#include <type_traits>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "include/initialization.h"
#include "include/lazy.h"
#include "include/naive.h"
#include "yinyang_simp.h" 
#include "quantized.h"
#include "nisk_kmeans.h"
#include "include/lsh.h"

template <typename T>
void update_centers_standard(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");


   parlay::parallel_for(0,k*d,[&] (size_t icoord) {

    size_t i = icoord/d; //the center we are looking at
    size_t coord = icoord % d; //coord we are looking at
    if (pts_grouped_by_center[i].second.size() > 0) {
      
      c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
      [&] (size_t ind) {
        return static_cast<float>(v[ind*d+coord]);
      }))/pts_grouped_by_center[i].second.size();

    }

   });

   t2.next("Added");



}

//nested parallel for loop improvement 
template <typename T>
void update_centers_2(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");


   parlay::parallel_for(0,k,[&] (size_t i) {

    parlay::parallel_for(0,d,[&] (size_t coord) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
        [&] (size_t ind) {
          return static_cast<float>(v[ind*d+coord]);
        }))/pts_grouped_by_center[i].second.size();

      }

    });

  
   });

   t2.next("Added");



}


//nested parallel for loop improvement 
template <typename T>
void update_centers_3(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");


   parlay::parallel_for(0,k,[&] (size_t i) {

    parlay::parallel_for(0,d,[&] (size_t coord) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
        [&] (size_t ind) {
          return static_cast<float>(v[ind*d+coord]);
        }))/pts_grouped_by_center[i].second.size();

      }

    },1); //does the 1 here help?

  
   });

   t2.next("Added");


}

//does the sort instead of a group by work faster?

// //nested parallel for loop improvement 
//unclear how to get the integer sorting to do what we want
// template <typename T>
// void update_centers_4(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
//   parlay::internal::timer t2;
//   t2.start();

//   parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

//   auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
//     return std::pair(asg[i],i);
//   }));

//   parlay::sequence<size_t> keys = parlay::tabulate(n,[&] (size_t i) {
//     return asg[i];
//   });

//   parlay::integer_sort_inplace(rangn,keys);
//   my_hist = parlay::histogram_by_key(keys);



//   t2.next("Grouped by");

//   abort();


//    parlay::parallel_for(0,k,[&] (size_t i) {

//     parlay::parallel_for(0,d,[&] (size_t coord) {
//        if (pts_grouped_by_center[i].second.size() > 0) {
        
//         c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
//         [&] (size_t ind) {
//           return static_cast<float>(v[ind*d+coord]);
//         }))/pts_grouped_by_center[i].second.size();

//       }

//     },1); //does the 1 here help?

  
//    });

//    t2.next("Added");



// }


//nested parallel for loop improvement 
//wow this much worse, 82s instead of 37s like the others
template <typename T>
void update_centers_5(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  std::cout << "method 5" << std::endl;
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");

    parlay::parallel_for(0,d,[&] (size_t coord) {
      parlay::parallel_for(0,k,[&] (size_t i) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
        [&] (size_t ind) {
          return static_cast<float>(v[ind*d+coord]);
        }))/pts_grouped_by_center[i].second.size();

      }

    }); 

  
   }); //does the 1 here help?

   t2.next("Added");


}


//nested parallel for loop improvement 
//also bad (79s for adding) <- appears that parallelizing by d first 
//makes something much worse (caching?)
template <typename T>
void update_centers_6(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  std::cout << "method 6" << std::endl;
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");

    parlay::parallel_for(0,d,[&] (size_t coord) {
      parlay::parallel_for(0,k,[&] (size_t i) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
        [&] (size_t ind) {
          return static_cast<float>(v[ind*d+coord]);
        }))/pts_grouped_by_center[i].second.size();

      }

    }); 

  
   }); //does the 1 here help?

   t2.next("Added");


}


//Can we improve the caching by pulling one point at a time?
template <typename T>
void update_centers_7(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  std::cout << "method 7" << std::endl;
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");

    parlay::parallel_for(0,k,[&] (size_t i) {
      size_t picked_center = pts_grouped_by_center[i].first;
      for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
        size_t point_coord = v[pts_grouped_by_center[i].second]
        for (size_t coord = 0; coord < d; coord++) {
          c[picked_center*d+coord] += static_cast<float>(v[pts_grouped_by_center[i].second])
        }
      }
      parlay::parallel_for(0,k,[&] (size_t i) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
        [&] (size_t ind) {
          return static_cast<float>(v[ind*d+coord]);
        }))/pts_grouped_by_center[i].second.size();

      }

    }); 

  
   }); //does the 1 here help?

   t2.next("Added");


}


//bench a given initialization method
//returns pair<double,double>: the time of the initialization method and the msse afterward
template<typename T, typename Initializer>         
void call_update_bench(T* v, size_t n, size_t d, size_t k, Distance& D, std::string output_folder) {


  float* c = new float[k*d];

  //Initializer needed before an update
  size_t* asg = new size_t[n];

  Initializer dummy; //just to get the name
  
  initialization_bench bencher(n,d,k,dummy.name());

  bencher.start_time();
  Initializer init;

  init(v,n,d,k,c,asg,D);

  bencher.end_time();  

  double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
  T* it2 = v+i*d;
  float buf[2048];
  for (unsigned short int j = 0; j < d; j++) {
  buf[j] = *it2;
  it2+=1;
  }
  return D.distance(buf,c+asg[i]*d,d);
  }))/n;//divide by n because mean sse

  bencher.set_msse(msse);

  //now onto the center benching

  parlay::internal::timer t2;
  t2.start();
  update_centers_6(v,n,d,k,c,asg);
  t2.next("Finished method"); //total time for completion

  delete[] c;
  delete[] asg;

  }



template <typename T>
void pick_init(T* v, size_t n, size_t d, size_t k, Distance& D, std::string init_choice, std::string output_folder) {
  if (init_choice == "None") {
    std::cout << "None picked aborting" << std::endl;
    abort();
  }
  else if (init_choice == "Standard") {
    //LSH guarantees each center has points
    call_update_bench<T,LazyStart<T>>(v,n,d,k,D,output_folder);
  }
 

}

int main(int argc, char* argv[]){
    commandLine P(argc, argv, "[-k <n_clusters>] [-i <input>] [-f <ft>] [-t <tp>] [-D <dist>]");

    size_t k = P.getOptionLongValue("-k", 10); // k is number of clusters
  
    std::string input = std::string(P.getOptionValue("-i", "")); // the input file
    std::string ft = std::string(P.getOptionValue("-f", "bin")); // file type, bin or vecs
    std::string tp = std::string(P.getOptionValue("-t", "uint8")); // data type
    std::string dist = std::string(P.getOptionValue("-D", "Euclidian")); // distance choice

    std::string init_choice = std::string(P.getOptionValue("-c", "None"));

    std::string output_folder = std::string(P.getOptionValue("-o",""));

    if(input == ""){ // if no input file given, quit
        std::cout << "Error: input file not specified" << std::endl;
        abort();
    }

    if((ft != "bin") && (ft != "vec")){ // if the file type chosen is not one of the two approved file types, quit 
    std::cout << "Error: file type not specified correctly, specify bin or vec" << std::endl;
    abort();
    }

    if((tp != "uint8") && (tp != "int8") && (tp != "float")){ // if the data type isn't one of the three approved data types, quit
        std::cout << "Error: vector type not specified correctly, specify int8, uint8, or float" << std::endl;
        abort();
    }

    if((ft == "vec") && (tp == "int8")){ // you can't store int8s in a vec file apparently I guess
        std::cout << "Error: incompatible file and vector types" << std::endl;
        abort();
    }

    // TODO: add support for vec files
    if (ft == "vec") {
        std::cout << "Error: vec file type not supported yet" << std::endl;
        abort();
    }

    Distance* D; // create a distance object, it can either by Euclidian or MIPS
    if (dist == "Euclidean") { 
        std::cout << "Using Euclidean distance" << std::endl;
        D = new EuclideanDistance();
    } else if (dist == "mips") {
        std::cout << "Using MIPS distance" << std::endl;
        D = new Mips_Distance();
    } else if (dist=="short") {
        std::cout << "Using short Euclidean" << std::endl;
        D = new EuclideanDistanceSmall();
    } 
    else {
        std::cout << "Error: distance type not specified correctly, specify Euclidean or mips" << std::endl;
        abort();
    }

    if (ft == "bin"){
        if (tp == "float") {
            auto [v, n, d] = parse_fbin(input.c_str());

            pick_init<float>(v,n,d,k,*D,init_choice,output_folder);
          
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            
            pick_init<uint8_t>(v,n,d,k,*D,init_choice,output_folder);
            
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            pick_init<int8_t>(v,n,d,k,*D,init_choice,output_folder);
            
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}

#endif //INIT_H