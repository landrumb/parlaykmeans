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



  void cluster(T* v, size_T n, size_t d, size_t k, float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
    size_t split = 4; //number of blocks to split dimensions into
    size_t len_per_split = d/4;
    size_t kstar = static_cast<size_t>(std::pow(k,1.0/split));
    if (split*(d/split) != d) {
      std::cout << "d not divisible by 4, aborting" << std::endl;
    }
    if (std::abs(std::pow(kstar,split)-k) >= .1) {
      std::cout << "k does not have 4th root, aborting" << std::endl;
    }
    //new point struct
    parlay::sequence<T*> vx(split);
    for (size_t i =0; i < vx.size(); i++) {
      vx[i] = new T[n*d/split];
    }

    parlay::parallel_for (0,n,[&] (size_t i) {
      for (size_t j = 0; j < d; j++) {
        vx[j/len_per_split][i*len_per_split+j%len_per_split] = v[i*d+j];
      }
    });
    parlay::sequence<kmeans_bench> loggers(split);

    parlay::sequence<float*> cx(split);
    for (size_t i = 0 ; i < split ; i++) {
      cx[i] = new float[kstar*len_per_split];
    }

    parlay::sequence<size_t*> asgx(split);
    for (size_t i = 0; i < split; i++) {
      asgx[i] = new size_t[n];
    }

    //new initialization, the previous initialization can't be copied?
    LazyStart<float> init;
    parlay::parallel_for(0,split,[&] (size_t i) {
      init(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D);
    });

    YinyangSimp<T> yy;

    parlay::parallel_for(0,split,[&] (size_t i) {
      loggers[i].start_time();
      yy.cluster(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D,loggers[i],max_iter,epsilon);
      loggers[i].end_time();
    });

    //copy answers into c, expanding outward
    for (size_t i = 0; i < k; i++) {
      parlay::sequence<size_t> center_selects(split);
      for (size_t j = 0; j < split; j++) {
        centers_selects[j] = (i / std::pow(kstar,j)) % std::pow(kstar,j+1);
      }
      for (size_t j = 0; j < d; j++) {
        c[i*d+j] = cx[j/len_per_split][center_selects[j/len_per_split]*d+j%len_per_split];
      }

    }

    //perhaps we ought to do one final assign step here
    //but for now just copying in whatever the 0th index has
    //alternatively could do a majority vote scheme
    //did this translation occur correctly?
    parallel_for(0,n,[&] (size_t i) {
      asg[i] = 0;
      for (size_t j = 0; j < split; j++) {
        asg[i] += asgx[j][i] * std::pow(kstar,j);
      }
    });




  }



};