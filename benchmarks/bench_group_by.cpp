//bench compare fast group by old group by, to see the improvement

#ifndef GROUP_BY_BENCH_H
#define GROUP_BY_BENCH_H

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


//helpful function for center calculation
//requires integer keys
template<typename index_type>
void fast_int_group_by(parlay::sequence<std::pair<index_type, parlay::sequence<index_type>>>& grouped, size_t n, index_type* asg) {
  parlay::internal::timer t3;
  t3.start();
   
   auto init_pairs = parlay::delayed_tabulate(n,[&] (size_t i) {
    return std::make_pair(asg[i],i);
  });

  t3.next("made pairs");
   parlay::sequence<std::pair<size_t,size_t>> int_sorted = parlay::integer_sort(init_pairs, [&] (std::pair<size_t,size_t> p) {return p.first;});
  //parlay::integer_sort_inplace(init_pairs, [&] (std::pair<size_t,size_t> p) {return p.first;});
   t3.next("int sorted");


    //store where each center starts in int_sorted
   auto start_pos = parlay::pack_index(parlay::delayed_tabulate(n,[&] (size_t i) {
    return i==0 || int_sorted[i].first != int_sorted[i-1].first;
  }));
  start_pos.push_back(n);

  t3.next("packed");
  grouped=parlay::tabulate(start_pos.size()-1, [&] (size_t i) {
    return std::make_pair(int_sorted[start_pos[i]].first, parlay::map(int_sorted.subseq(start_pos[i],start_pos[i+1]),[&] (std::pair<size_t,size_t> ind) { return ind.second;}) );

  });
  t3.next("grouped");



}

//bench a given initialization method
//returns pair<double,double>: the time of the initialization method and the msse afterward
template<typename T, typename Initializer>         
void bench_assign(T* v, size_t n, size_t d, size_t k, Distance& D, std::string output_folder) {

  float* c = new float[k*d];
  size_t* asg = new size_t[n];

  Initializer dummy; //just to get the name
  
  initialization_bench bencher(n,d,k,dummy.name());

  bencher.start_time();
  Initializer init;

  init(v,n,d,k,c,asg,D);

  bencher.end_time();  

  //onto group_by

  parlay::internal::timer t2;
  t2.start();
  auto rangn = parlay::delayed_tabulate(n,[&] (size_t i) {
    return i;
  });
  parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
    }));

  double standard_group_time = t2.next_time();

  parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center2;
  fast_int_group_by(pts_grouped_by_center2,n,asg);
  

  
  double float_assign_time = t2.next_time();

  std::cout << "group by time: standard, fast \t" << standard_group_time << " " << float_assign_time << std::endl;

 // std::cout << "assign times, T and float resp.: " << t_assign_time << ", " << float_assign_time << std::endl;
  

  delete[] c;
  delete[] asg;

  // std::cout << init.name() << ": time-- " << init_time << ", msse-- " << msse << std::endl;
  // std::cout << "n: " << n << ", d: " << d << ", k: " << k << ", Distance : " << D.id() << std::endl;
  //return std::make_pair(init_time,msse);

  }

template <typename T>
void pick_assign(T* v, size_t n, size_t d, size_t k, Distance& D, std::string init_choice, std::string output_folder) {
  if (init_choice == "None") {
    std::cout << "None picked aborting" << std::endl;
    abort();
  }
  else if (init_choice == "Standard") {
    bench_assign<T,LazyStart<T>>(v,n,d,k,D,output_folder);
  }

}

int main(int argc, char* argv[]){
    commandLine P(argc, argv, "[-k <n_clusters>] [-i <input>] [-f <ft>] [-t <tp>] [-D <dist>] [-n <chosen_k>]");

    size_t k = P.getOptionLongValue("-k", 10); // k is number of clusters
    size_t chosen_n = P.getOptionLongValue("-n",1000);

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

            pick_assign<float>(v,chosen_n,d,k,*D,init_choice,output_folder);
          
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            
            pick_assign<uint8_t>(v,chosen_n,d,k,*D,init_choice,output_folder);
            
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            pick_assign<int8_t>(v,chosen_n,d,k,*D,init_choice,output_folder);
            
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}

#endif //ASSIGN_H