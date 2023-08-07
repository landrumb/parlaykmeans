#ifndef LSHQ_H
#define LSHQ_H

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

#include "yinyang_simp.h" //can switch to fast_center
#include "lsh.h"


//the other quantized.h file is doing a cartesian k-means approach
//whereas this version is using PQ (via LSH) to reduce the number
//of dimensions we are looking at

template<typename T>
struct LSHQuantizedKmeans {

  // //TODO improve the bit set distance calculation!
  // for later?
  // uint8_t bitset_dist(const std::bitset& a, const std::bitset& b, const size_t& num_hp) {
  //   uint8_t dist_sum = 0;
  //   for (uint8_t i = 0; i < num_hp; i++) {
  //     if (a[i] != b[i]) dist_sum += 1;
  //   }

  // }

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {

    parlay::internal::timer t2;
    t2.start();

    //the number of hyperplanes we are looking at must be BITSET_MAX in the LSH file for the LSH code to work
    const size_t num_hp = 256;
    size_t new_d = 32;

    LSH<T> lsh_runner;

    parlay::sequence<parlay::sequence<float>> hps(num_hp,parlay::sequence<float>(d));

    lsh_runner.generate_hps(hps,d,num_hp);


    uint8_t* pqv = new uint8_t[n*new_d];
    float* pqc_float = new float[k*new_d];



   parlay::parallel_for(0,n,[&] (size_t i) {
    lsh_runner.template get_uint8_reduced_pt<T,uint8_t>(v+i*d,pqv+i*new_d,hps,d,num_hp);

   });

  //TODO can we use the same function to product quant our initial centers? 
  //knowing they'll be cast to floats
  //think so?
   parlay::parallel_for(0,k,[&] (size_t i) {
    lsh_runner.template get_uint8_reduced_pt<float,float>(c+i*d,pqc_float+i*new_d,hps,d,num_hp);
   });

   std::cout << "Debugging, printing new centers "<< std::endl;
   for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      std::cout << pqc_float[i*new_d+j] << " ";
    }
    std::cout << std::endl;
   }

  //concerned about distance function for d=36 or smaller
  Distance* D2 = new EuclideanDistanceSmall();

  kmeans_bench other_logger = kmeans_bench(n,new_d,k,max_iter, epsilon,"Lazy","Internal YY to LSHQ");

  //use yy as a subroutine
  YinyangSimp<uint8_t> yy;
   yy.cluster(pqv,n,new_d,k,pqc_float,asg,*D2,logger,max_iter,epsilon);

   //from the assignments, build the centers
   //

   parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t i) {
    return i;
   });

  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
      return std::pair(asg[i],i);
    }));

   parlay::parallel_for(0,k*d,[&] (size_t icoord) {

    size_t i = icoord/d; //the center we are looking at
    size_t coord = icoord % d; //coord we are looking at
    if (pts_grouped_by_center[i].second.size() > 0) {
      //TODO FIXME .second[ind] is doubling ind, should just be ind here??
      //because a map takes out the actual values from the sequence
      c[pts_grouped_by_center[i].first*d+coord] = parlay::reduce(parlay::map(pts_grouped_by_center[i].second,
      [&] (size_t ind) {
        return static_cast<float>(v[pts_grouped_by_center[i].second[ind]*d+coord]);
      }))/pts_grouped_by_center[i].second.size();

    }

   });
    double run_time = t2.next_time();

   double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
    T* it2 = v+i*d;
    float buf[2048];
    for (unsigned short int j = 0; j < d; j++) {
      buf[j] = *it2;
      it2+=1;
    }
    return D.distance(buf,c+asg[i]*d,d);
  }));

      parlay::sequence<float> zero_filled_seq(k,0);

   logger.add_iteration(run_time,0,msse,0,0,zero_filled_seq,t2.next_time());

  }
};

#endif //LSHQ_H