//use LSH bucketing for a good initialization
//unclear if this has been done before
//TODO: check
#ifndef LSHCODE
#define LSHCODE

#include "include/utils/NSGDist.h"

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
#include <bitset>


template <typename T>
struct LSH {

  static const size_t BITSET_MAX = 32;

  std::string name() {
    return "LSH INIT";
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

  void operator()(T* v, size_t n,size_t d,size_t k, float* c, size_t* asg, Distance& D) {

    std::cout << "starting lsh" << std::endl;
     parlay::sequence<parlay::sequence<float>> hps(BITSET_MAX,parlay::sequence<float>(d));
     generate_hps(hps,d);
     
     //<hash value, point id> pair
     parlay::sequence<std::pair<size_t, size_t>> pts_with_hash = parlay::tabulate(n, [&] (size_t i) {
      return std::make_pair(get_hash(v+i*d,hps,n,d),i);
     });

         std::cout << "boutaa sort" << std::endl;

      //TODO how to use a parlay sort? 
      //getting the below error
      // error: expected primary-expression before ')' token
   //63 |      parlay::sort_inplace(pts_with_hash,less_pair);
   //oh need to pass an object not a type**
      parlay::sort_inplace(pts_with_hash,less_pair());
      //std::sort(pts_with_hash.begin(),pts_with_hash.end(),less_pair());

     if (n%k != 0) {
      std::cout << "uneven bucketing, aborting " << std::endl;
      abort();
     }

       std::cout << "calculating centers" << std::endl;

     //after sorting by LSH hash we can assign
     //yikes taking a while not well parallelized 
    //  parlay::parallel_for(0,k, [&] (size_t i) {

    //   size_t id = i*d;
    //   //noverk points per bucket
    //   size_t noverk = n/k;

    //   //make sure c is set to 0
    //   for (size_t j = id; j < (i+1)*d; j++) {
    //     c[j] = 0;
    //   }

      

    //   for (size_t j = i * noverk; j <  (i+1) * noverk; i++) {
    //     asg[pts_with_hash[j].second] = i;
    //     for (size_t coord = 0; coord < d; coord++) {
    //       c[id+coord] += v[pts_with_hash[j].second*d+coord];
    //     }
    //   }

    //   for (size_t coord = 0; coord < d; coord++) {
    //     c[id+coord] /= noverk;
    //   }


    //  },1); //granularity 1 needed? TODO CHECK

    //try again:
    //ohhh use a reduce!
    //woah the reduce much faster
    size_t noverk = n/k;

    parlay::sequence<size_t> rang_noverk = parlay::tabulate(noverk,[&] (size_t i) {
      return i;

    });

     parlay::parallel_for(0,k*d, [&] (size_t icoord) {

      //TODO remove repeated division, this is extra work
      size_t i = icoord / d; //which center we look at
      size_t coord = icoord % d; //which coord we are accessing

      //make sure c is set to 0
      c[icoord] = parlay::reduce(parlay::map(rang_noverk,[&] (size_t j) {
        return static_cast<float>(v[pts_with_hash[j+i*noverk].second*d+coord]);
      }))/noverk;
      
     });

     std::cout << "Finished with center added, now onto asg" << std::endl;

    //bad way to do it
    //  parlay::parallel_for(0,k, [&] (size_t i) {
    //   //noverk points per bucket
    //   size_t noverk = n/k;

    //   for (size_t j = i * noverk; j <  (i+1) * noverk; i++) {
    //     asg[pts_with_hash[j].second] = i;
        
    //   }

    //  }); 

     parlay::parallel_for(0,n,[&] (size_t i) {
      size_t chosen_center = i / noverk;
      asg[pts_with_hash[i].second] = chosen_center;

     });

     std::cout << "left lsh code " << std::endl;

  }

  void generate_hps(parlay::sequence<parlay::sequence<float>>& hps, size_t d) {
    std::random_device r;
    std::default_random_engine rng{r()};
    std::normal_distribution<float> gaussian_dist;
    for (size_t i = 0; i < BITSET_MAX; i++) {
      hps[i] = parlay::sequence<float>(d);
      for (size_t j = 0; j < d; j++) {
        hps[i][j] = gaussian_dist(rng);
      }
    }

  }

  //given a point and the hyperplanes, get the hash id of the pt
  size_t get_hash(T* pt, const parlay::sequence<parlay::sequence<float>>& hps, size_t n, size_t d) {
    std::bitset<BITSET_MAX> bin_list; 
    //rang gives 0-31
    parlay::sequence<size_t> rang = parlay::tabulate(BITSET_MAX, [&] (size_t i) {return i;});

    for (size_t i = 0; i < BITSET_MAX; i++) {
      //TODO does the float cast on pt happen automatically? if so remove
      //the manual cast
      float dot_p = parlay::reduce(parlay::map(rang, [&] (size_t j) {return hps[i][j] * static_cast<float>(pt[j]);}));
      if (dot_p > 0) {
        bin_list[i] = 1;
      }
      else {
        bin_list[i] = 0;
      }
    }
    //to_ulong a bitset method to get a long out of a bitset
    return bin_list.to_ulong();

  }

};

#endif //LSHCODE