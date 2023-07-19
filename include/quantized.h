#include <cmath>
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
#include "include/utils/kmeans_bench.h"
#include "yinyang_simp.h"
//#include "include/utils/union_find.h"

//use quantized version to reduce runtimes
template <typename T>
struct QuantizedKmeans {



  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {

    parlay::internal::timer t2;
    t2.start();

    size_t split = 4; //number of blocks to split dimensions into
    size_t len_per_split = d/split;
    size_t kstar = static_cast<size_t>(std::pow(k,1.0/split));
    if (split*len_per_split != d) {
      std::cout << "d not divisible by 4, aborting" << std::endl;
    }
    if (std::abs(std::pow(kstar,split)-k) >= .1) {
      std::cout << "k does not have 4th root, aborting" << std::endl;
    }
    //new point struct
    //TODO make kmeans method that doesn't need this copy
    parlay::sequence<T*> vx(split);
    for (size_t i =0; i < vx.size(); i++) {
      vx[i] = new T[n*d/split];
    }

    parlay::parallel_for (0,n,[&] (size_t i) {
      for (size_t j = 0; j < d; j++) {
        vx[j/len_per_split][i*len_per_split+j%len_per_split] = v[i*d+j];
      }
    });

    //Debug! need to initialize? () init not allowed
    parlay::sequence<kmeans_bench> loggers;
    for (size_t i = 0 ; i < split; i++) {
      loggers.push_back(kmeans_bench(n,len_per_split,kstar,max_iter,epsilon,"Lazy","Piece of Quantized"));
    }

    parlay::sequence<float*> cx(split);
    for (size_t i = 0 ; i < split ; i++) {
      cx[i] = new float[kstar*len_per_split];
    }

    parlay::sequence<size_t*> asgx(split);
    for (size_t i = 0; i < split; i++) {
      asgx[i] = new size_t[n];
    }

    parlay::sequence<float> zero_filled_seq(0);
    logger.add_iteration(0,0,51,0,0,zero_filled_seq,t2.next_time());

    //new initialization, the previous initialization can't be copied?
    //TODO how to use previous init?
    LazyStart<T> init;
    parlay::parallel_for(0,split,[&] (size_t i) {
      init(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D);
    });

    logger.add_iteration(0,0,52,0,0,zero_filled_seq,t2.next_time());

    YinyangSimp<T> yy;

    //Because printing in parallel a mess, suppresses internal yy logging
    parlay::parallel_for(0,split,[&] (size_t i) {
      loggers[i].start_time();
      yy.cluster(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D,loggers[i],max_iter,epsilon,true);
      loggers[i].end_time();
    });

    logger.add_iteration(t2.next_time(),0,53,0,0,zero_filled_seq,0);


    //copy answers into c, expanding outward
    for (size_t i = 0; i < k; i++) {
      parlay::sequence<size_t> center_selects(split);
      //std::cout << "printing center select for i: " << i << std::endl;

      for (size_t j = 0; j < split; j++) {
        center_selects[j] = static_cast<long>(i / std::pow(kstar,j)) % kstar;
        //std::cout << center_selects[j] << " ";
      }
      //std::cout << std::endl;
    
      for (size_t j = 0; j < d; j++) {
        c[i*d+j] = cx[j/len_per_split][center_selects[j/len_per_split]*d+j%len_per_split];
      }

    }

    //perhaps we ought to do one final assign step here
    //but for now just copying in whatever the 0th index has
    //alternatively could do a majority vote scheme
    //did this translation occur correctly?
    parlay::parallel_for(0,n,[&] (size_t i) {
      asg[i] = 0;
      for (size_t j = 0; j < split; j++) {
        asg[i] += asgx[j][i] * std::pow(kstar,j);
      }
    });

    //Check msse
    double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
      return D.distance(v+i*d,c+asg[i]*d,d);
    }));
    msse /= n; //divide by n because mean msse


    logger.add_iteration(0,t2.next_time(),msse,0,0,zero_filled_seq,0);

  }



};