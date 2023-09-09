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
//nice!! this one does the add in just 3 seconds by doing the caching right
//for k=100K works nice, I wonder if this sucks for small k (k=10?)
//Wow k=10 it's much worse (97s)! 
//how does scheduling work such that higher k reduces time significantly? 
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
      size_t picked_center_d = pts_grouped_by_center[i].first*d;
      for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
        size_t point_coord = pts_grouped_by_center[i].second[j]*d;
        for (size_t coord = 0; coord < d; coord++) {
          c[picked_center_d + coord] += static_cast<float>(v[point_coord + coord]);
        }
      }
      

  });

  t2.next("Added");


   parlay::parallel_for(0,k,[&] (size_t i) {

    parlay::parallel_for(0,d,[&] (size_t coord) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] /= pts_grouped_by_center[i].second.size();

      }

    });
  
   });

   t2.next("Divided");


}


  //comparator for pairs of size_t
  //compare the first elt then the second elt
  struct less_pair {
    bool operator() (const std::pair<size_t,size_t>& x, const std::pair<size_t,size_t>& y) const {
      if (x.first < y.first) {
        return true;
      }
      else if (x.first == y.first) {
        return x.second < y.second;
      }
      else {
        return false;
      }
      
    }

     bool comp(const std::pair<size_t,size_t>& x, const std::pair<size_t,size_t>& y) const {
      if (x.first < y.first) {
        return true;
      }
      else if (x.first == y.first) {
        return x.second < y.second;
      }
      else {
        return false;
      }
      
    }
  };

template <typename T>
void my_group_by(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg,
parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& grouped) {

  parlay::internal::timer t3;
  t3.start();


  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  parlay::sequence<std::pair<size_t,size_t>> paired = parlay::tabulate(n,[&] (size_t i) {
    return std::make_pair(asg[i],i);
  });

  t3.next("paired");
  
  
  parlay::sort_inplace(paired,less_pair());

  t3.next("Sorted");

  parlay::sequence<size_t> changed = parlay::tabulate(n,[&] (size_t i) {
    if (i==0) return static_cast<size_t>(0);
    else {
      if (paired[i].first != paired[i-1].first) {
        return static_cast<size_t>(1);
      }
      else {
        return static_cast<size_t>(0);
      }

    }
  });

  t3.next("changed");
  parlay::sequence<size_t> points_of_change = parlay::filter(rangn,[&] (size_t i) {
    return (i==0 || changed[i] == 1);
  });
    t3.next("filtered");

  auto [offsets, total] = parlay::scan(points_of_change);

  t3.next("scanned");
  std::cout << "about to allo" << std::endl;

  parlay::parallel_for(0,k,[&] (size_t i) {
    if (i==k-1) {
      grouped[i].second = parlay::sequence<size_t>(n-offsets[i]);
    }
    else {
      grouped[i].second = parlay::sequence<size_t>(offsets[i+1]-offsets[i]);
    }
  });

  t3.next("allocated space");


  parlay::parallel_for(0,k,[&] (size_t i) {
    if (i == k-1) {
      for (size_t j = offsets[i]; j < n; j++) {
        grouped[i].second.push_back(paired[j].second);
      }
    }
    else {
      for (size_t j = offsets[i]; j < offsets[i+1]; j++) {
        grouped[i].second.push_back(paired[i].second);


      }

    }
  
    
  });

  t3.next("Pushed back");

 

}


template <typename T>
void update_centers_8(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  std::cout << "method 8" << std::endl;
  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));

  t2.next("Grouped by");



  parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> grouped = parlay::tabulate(k,[&] (size_t i) {
    return std::make_pair(i,parlay::sequence<size_t>());
  });

  my_group_by(v,n,d,k,c,asg,grouped);

  t2.next("My grouped by");


  abort(); //just case about group by

    parlay::parallel_for(0,k,[&] (size_t i) {
      size_t picked_center_d = pts_grouped_by_center[i].first*d;
      for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
        size_t point_coord = pts_grouped_by_center[i].second[j]*d;
        for (size_t coord = 0; coord < d; coord++) {
          c[picked_center_d + coord] += static_cast<float>(v[point_coord + coord]);
        }
      }
      

  });

  t2.next("Added");


   parlay::parallel_for(0,k,[&] (size_t i) {

    parlay::parallel_for(0,d,[&] (size_t coord) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] /= pts_grouped_by_center[i].second.size();

      }

    });
  
   });

   t2.next("Divided");


}

template <typename T>
void get_groups1(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& pts_grouped_by_center) {


  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) { return i; });
  pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
  }));
}

template <typename T>
void add_points1(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& pts_grouped_by_center) {

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


}

template <typename T>
void add_points2(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>>& pts_grouped_by_center) {

  parlay::parallel_for(0,k,[&] (size_t i) {
      size_t picked_center_d = pts_grouped_by_center[i].first*d;
      for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
        size_t point_coord = pts_grouped_by_center[i].second[j]*d;
        for (size_t coord = 0; coord < d; coord++) {
          c[picked_center_d + coord] += static_cast<float>(v[point_coord + coord]);
        }
      }
      

  },1); //does gran 1 help here? 

    parlay::parallel_for(0,k,[&] (size_t i) {

    parlay::parallel_for(0,d,[&] (size_t coord) {
       if (pts_grouped_by_center[i].second.size() > 0) {
        
        c[pts_grouped_by_center[i].first*d+coord] /= pts_grouped_by_center[i].second.size();

      }

    });
  
   });

}

template<typename T>
std::pair<double,double> update_centers_segmented(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg) {
  parlay::internal::timer t2;
  t2.start();
  parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center;
  get_groups1(v,n,d,k,c,asg,pts_grouped_by_center);
  double group_time = t2.next_time();
  add_points1(v,n,d,k,c,asg,pts_grouped_by_center);
  double add_time = t2.next_time();
  return std::make_pair(group_time,add_time);
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
  update_centers_8(v,n,d,k,c,asg);
  t2.next("Finished method"); //total time for completion

  delete[] c;
  delete[] asg;

  }

//benches just the adding phase
template<typename T, typename Initializer>
std::tuple<double,double,double> call_adding_bench(T* v, size_t n, size_t d, size_t k, Distance& D, std::string output_folder) {


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

  parlay::internal::timer t2;
  t2.start();

  parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center;
  get_groups1(v,n,d,k,c,asg,pts_grouped_by_center);

  double group_time = t2.next_time();
  std::cout << "Grouped in "<< group_time << std::endl;

  //now onto the center benching
  //just the add step

  add_points1(v,n,d,k,c,asg,pts_grouped_by_center);
  double add1 = t2.next_time();
  std::cout << "method 1 for adding points " << add1 << std::endl;
  add_points2(v,n,d,k,c,asg,pts_grouped_by_center);
  double add2 = t2.next_time();
  std::cout << "method 2 for adding points " << add2 << std::endl;

 
  // std::pair<double, double> time_output = update_centers_segmented(v,n,d,k,c,asg);
  
  // std::cout << "Group time: " << time_output.first << ", Add time: " << time_output.second << std::endl;

  delete[] c;
  delete[] asg;

  return std::make_tuple(group_time,add1,add2);

}


//benches just the adding phase
//over n =1mil, 10mil, 100mil, 200mil,400mil,600mil,800mil,1bil
template<typename T, typename Initializer>
void call_adding_bench_overn(T* v, size_t n, size_t d, size_t k, Distance& D, std::string output_folder) {

  parlay::sequence<size_t> values_of_n = {1'000'000,10'000'000,100'000'000, 200'000'000,400'000'000,600'000'000,800'000'000,1'000'000'000};
  parlay::sequence<std::tuple<double,double,double>> results(values_of_n.size());

  //note we're doing the runs sequentially so they don't step on each other!
  for (size_t i = 0; i < values_of_n.size(); i++) {
    results[i] = call_adding_bench<T,Initializer>(v,values_of_n[i],d,k,D,"");

  }

  //output to csv
  if (output_folder != "") {
     auto c_t = std::time(nullptr);
  auto tm = *std::localtime(&c_t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%m_%d_%Y");
  auto date_str = oss.str();
  std::string fname = output_folder + "data.bench_update_add_" + std::to_string(k) + "_" + date_str + ".csv";
  std::cout << "outputting to csv: " << fname << std::endl;
  std::ofstream file(fname);
  file << "d" << ", " << d << std::endl;
  file << "k" << ", " << k << std::endl;
  
  for (size_t i = 0; i < values_of_n.size(); i++) {
    file << "n" << ", " << values_of_n[i] << ", " << std::get<0>(results[i]) << ", " << std::get<1>(results[i]) << ", " << std::get<2>(results[i]) << std::endl;

  }


  }
 

}



template<typename T, typename Initializer>
void pick_init_stage2(T* v, size_t n, size_t d, size_t k, Distance& D, std::string output_folder, std::string bench_mode) {
  if (bench_mode == "normal") {
    call_update_bench<T,Initializer>(v,n,d,k,D,output_folder);
  }
  else if (bench_mode == "add_only") {
    call_adding_bench<T,Initializer>(v,n,d,k,D,output_folder);


  }
  else if (bench_mode == "add_overn") {
    call_adding_bench_overn<T,Initializer>(v,n,d,k,D,output_folder);
  }
  else {
    std::cout << "Did not pick bench mode aborting" << std::endl;
    abort();
  }

}

//separating into nested i/o functions to reduce if size
template <typename T>
void pick_init(T* v, size_t n, size_t d, size_t k, Distance& D, std::string init_choice, std::string output_folder, std::string bench_mode) {
  if (init_choice == "None") {
    std::cout << "None picked aborting" << std::endl;
    abort();
  }
  else if (init_choice == "Standard") {
    //LSH guarantees each center has points
    pick_init_stage2<T,LazyStart<T>>(v,n,d,k,D,output_folder,bench_mode);
  }
 
}


int main(int argc, char* argv[]){
    commandLine P(argc, argv, "[-k <n_clusters>] [-i <input>] [-f <ft>] [-t <tp>] [-D <dist>] [-n <n_points>] [-v <which bench type>]");

    size_t k = P.getOptionLongValue("-k", 10); // k is number of clusters
  
    std::string input = std::string(P.getOptionValue("-i", "")); // the input file
    std::string ft = std::string(P.getOptionValue("-f", "bin")); // file type, bin or vecs
    std::string tp = std::string(P.getOptionValue("-t", "uint8")); // data type
    std::string dist = std::string(P.getOptionValue("-D", "Euclidian")); // distance choice

    std::string init_choice = std::string(P.getOptionValue("-c", "None"));

    std::string output_folder = std::string(P.getOptionValue("-o",""));

    std::string bench_mode = std::string(P.getOptionValue("-v",""));

    size_t truen = P.getOptionLongValue("-n",1000000);

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


            pick_init<float>(v,truen,d,k,*D,init_choice,output_folder, bench_mode);
            
          
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            
            pick_init<uint8_t>(v,truen,d,k,*D,init_choice,output_folder,bench_mode);
            
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            pick_init<int8_t>(v,truen,d,k,*D,init_choice,output_folder,bench_mode);
            
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}

#endif //INIT_H