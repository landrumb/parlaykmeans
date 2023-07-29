//file to rate the quality of different initialization methods

#ifndef INIT_H
#define INIT_H

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

#include "include/initialization.h"
#include "include/lazy.h"
#include "include/naive.h"
#include "yinyang_simp.h" 
#include "quantized.h"
#include "nisk_kmeans.h"
#include "include/lsh.h"

//bench a given initialization method
//returns pair<double,double>: the time of the initialization method and the msse afterward
template<typename T, typename Initializer>         
std::pair<double,double> bench_init(T* v, size_t n, size_t d, size_t k, Distance& D) {

  parlay::internal::timer t;

  float* c = new float[k*d];
  size_t* asg = new size_t[n];

  std::cout << "starting bench" << std::endl;
  
  t.start();

  Initializer init;
  init(v,n,d,k,c,asg,D);

  double init_time = t.next_time();

  

  //Check msse
  double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
    T* it2 = v+i*d;
    float buf[2048];
    for (unsigned short int j = 0; j < d; j++) {
      buf[j] = *it2;
      it2+=1;
    }
    return D.distance(buf,c+asg[i]*d,d);
  }));
  msse /= n; //divide by n because mean sse

  delete[] c;
  delete[] asg;

  std::cout << init.name() << ": time-- " << init_time << ", msse-- " << msse << std::endl;
  std::cout << "n: " << n << ", d: " << d << ", k: " << k << ", Distance : " << D.id() << std::endl;
  return std::make_pair(init_time,msse);
}
template <typename T>
void pick_init(T* v, size_t n, size_t d, size_t k, Distance& D, std::string init_choice) {
  if (init_choice == "None") {
    std::cout << "None picked aborting" << std::endl;
    abort();
  }
  else if (init_choice == "Lazy") {
    bench_init<T,LazyStart<T>>(v,n,d,k,D);
  }
  else if (init_choice == "MacQueen") {
    bench_init<T,MacQueen<T>>(v,n,d,k,D);
  }
  else if (init_choice == "Forgy") {
    bench_init<T,Forgy<T>>(v,n,d,k,D);

  }
  else if (init_choice == "KmeansPlusPlus") {
    bench_init<T,KmeansPlusPlus<T>>(v,n,d,k,D);
  }
  else if (init_choice == "LSH") {
    bench_init<T,LSH<T>>(v,n,d,k,D);
  }
  else {
    std::cout << "aborting, wrong spelling " << init_choice << std::endl;
    abort();
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

            pick_init<float>(v,n,d,k,*D,init_choice);
          
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            
            pick_init<uint8_t>(v,n,d,k,*D,init_choice);
            
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            pick_init<int8_t>(v,n,d,k,*D,init_choice);
            
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}

#endif //INIT_H