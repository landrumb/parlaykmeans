#ifndef QUANT
#define QUANT

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
#include "include/nisk_kmeans.h"
//#include "include/utils/union_find.h"

//use quantized version to reduce runtimes
template <typename T>
struct QuantizedKmeans {

  //for my sanity
  typedef parlay::sequence<parlay::sequence<float>> matrix;
  typedef parlay::sequence<parlay::sequence<T>> tmat;


  parlay::sequence<size_t> rand_seq(size_t n, size_t k){
      assert(n > k);
      std::random_device randomizer;
      std::mt19937 generate(randomizer());
    // std::uniform_int_distribution<size_t> dis(1, n);

      parlay::sequence<size_t> random_numbers(k);
      size_t bucket = n / k;

      parallel_for(0, k, [&](size_t i){
          std::uniform_int_distribution<size_t> dis(bucket * i, bucket * (i+1) - 1);
          random_numbers[i] = dis(generate);
      });
    

    return random_numbers;
  }


  //given a data sequence and a codebook, return an encoded version
  parlay::sequence<size_t> quantize_pq(parlay::sequence<T>& data, 
  parlay::sequence<parlay::sequence<parlay::sequence<float>>>& subcenters,
  size_t d, size_t m, size_t h, Distance& D) {
    parlay::sequence<size_t> code(m,0);
    //copy point into buffer
    float buf[2048];
    T* it = parlay::make_slice(data).begin();
    for (size_t j = 0; j < d; j++) {
        buf[j]=*it;
        it += 1; //increment the pointer
    }
    //for each block, add the id of our closest subcenter to our code
    for (size_t i = 0; i < m; i++) {
      auto distances = parlay::delayed::map(subcenters[i], [&] (parlay::sequence<float> c) {
        return D.distance(parlay::make_slice(c).begin(),buf,d);
      });
      //is this parlay::min_elt? TODO where is min_elt coming from, std?
      code[i] = std::min_element(distances) - distances.begin();
    }
    return code;
  }

  //quantization method will add here
  //data <- data vector we want to quantize
  //R <- rotation matrix 
  //d <- # of dim
  //m <- # of blocks of subcenters
  //h <- # of subcenters per block
  parlay::sequence<size_t> quantize_opq(parlay::sequence<T>& data, 
  parlay::sequence<parlay::sequence<float>>& R, 
  parlay::sequence<parlay::sequence<parlay::sequence<float>>>& subcenters,
  size_t d, size_t m, size_t h, Distance& D) {
    parlay::sequence<T> new_data = parlay::tabulate(d,[&] (size_t i) {
      return parlay::reduce(parlay::tabulate(d,[&] (size_t j) {
        return data[j] * R[i][j];
      }));
    });
    return quantize_pq(new_data,subcenters,d,m,h,D);
  }

  //train the product quant
  //returns:  
  //subcenters <- list of the subcenter matrices for each block and
  //R (rotation matrix)
  //using an identity (natural) initialization
  std::pair<matrix, matrix> train_opq (tmat pts, size_t n, size_t d, size_t m, size_t h, size_t niter) {

    if (d % m != 0) {
      std::cout << "does not split evenly d, m, aborting\n";
      abort();
    }
    size_t subdim = d/m;

    matrix R(d,parlay::sequence<float>(d));
    for (size_t i =0 ; i  < d; i++) {
      R[i][i]=1;
    }
    parlay::sequence<matrix> subcenters(m,matrix(d,parlay::sequence<float>(h)));

    for (size_t i = 0; i < m; i++) {
      parlay::sequence<size_t> perm = rand_seq(n,h);
      //grabbing a reference for ease of typing
      matrix& current_center = subcenters[i];
      for (size_t j = 0; j < h; j++) {
        size_t chsn_pt_coord = perm[j]; //chosen point coordinate
        //note that centers go down columns
        for (size_t coord = i*subdim; coord < (i+1)*subdim; coord++) {
          current_center[coord][j] = pts[coord][chsn_pt_coord];

        }

      }
    }

    //do clustering

    //do svd to get out rotation matrix
    //how do we do svd? is it built into parlay? use eigen or something like that?

    
  }



  
  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {

    parlay::internal::timer t2;
    t2.start();

    size_t split = 4; //number of blocks to split dimensions into
    size_t len_per_split = d/split;
    size_t kstar = static_cast<size_t>(std::pow(k,1.0/split));

    std::cout << "kstar is " << kstar << std::endl;
    std::cout << "len_per_split is " << len_per_split << std::endl;

    if (split*len_per_split != d) {
      std::cout << "d not divisible by 4, aborting" << std::endl;
      abort();

    }
    if (std::abs(std::pow(kstar,split)-k) >= .1) {
      std::cout << "k does not have 4th root, aborting" << std::endl;
      abort();
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
      loggers.push_back(kmeans_bench(n,len_per_split,kstar,max_iter,epsilon,"Lazy","Piece of Quantized YY"));
    }

    parlay::sequence<float*> cx(split);
    for (size_t i = 0 ; i < split ; i++) {
      cx[i] = new float[kstar*len_per_split];
    }

    parlay::sequence<size_t*> asgx(split);
    for (size_t i = 0; i < split; i++) {
      asgx[i] = new size_t[n];
    }

    parlay::sequence<float> zero_filled_seq(k,0);
    logger.add_iteration(0,0,51,0,0,zero_filled_seq,t2.next_time());

    //new initialization, the previous initialization can't be copied?
    //TODO how to use previous init?
    LazyStart<T> init;
    parlay::parallel_for(0,split,[&] (size_t i) {
    //for (size_t i = 0; i < split; i++) {
      init(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D);
    //}
    });

    logger.add_iteration(0,0,52,0,0,zero_filled_seq,t2.next_time());

    NiskKmeans<T> nisk;

    //Because printing in parallel a mess, suppresses internal yy logging
    parlay::parallel_for(0,split,[&] (size_t i) {
   // for (size_t i = 0; i < split ; i++) {
    
      loggers[i].start_time();
      nisk.cluster(vx[i],n,len_per_split,kstar,cx[i],asgx[i],D,loggers[i],max_iter,epsilon);
      loggers[i].end_time();
      std::cout << "Finished a split" << std::endl;
   // }
   });

    std::cout << "Finished split kmeans" << std::endl;

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

    logger.add_iteration(0,t2.next_time(),54,0,0,zero_filled_seq,0);


    //perhaps we ought to do one final assign step here
    //but for now just copying in whatever the 0th index has
    //alternatively could do a majority vote scheme
    //did this translation occur correctly?
    //TODO FIXME
    parlay::parallel_for(0,n,[&] (size_t i) {
      asg[i] = 0;
      for (size_t j = 0; j < split; j++) {
        asg[i] += asgx[j][i] * std::pow(kstar,j);
      }
    });


    //Check msse
    //Whoops need to copy point to float for calc?
    double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
      T* it2 = v+i*d;
      float buf[2048];
      for (unsigned short int j = 0; j < d; j++) {
        buf[j] = *it2;
        it2+=1;
      }
      return D.distance(buf,c+asg[i]*d,d);
    }));
    msse /= n; //divide by n because mean msse

    //This msse calculate version takes a really long time
    // //rang is range from [0,k)
    // parlay::sequence<size_t> rang = parlay::tabulate(k,[&] (size_t i) {
    //   return i;
    // });

    // //calculate msse and do assignments together
    // double msse = parlay::reduce(parlay::tabulate(n,[&] (size_t i) {
    //   float buf[2048];
    //   T* it = v+i*d;
    //   for (size_t j = 0; j < d; j++) {
    //       buf[j]=*it;
    //       it += 1; //increment the pointer
    //   }
    //   auto distances = parlay::delayed::map(rang, [&] (size_t k) {
    //     return D.distance(buf,c+k*d,d);
    //   });
    //   asg[i] = parlay::min_element(distances) - distances.begin();
    //   if (asg[i] >= k || asg[i] < 0) {
    //     std::cout << "seg fault leave" << std::endl;
    //   }
    //   return distances[asg[i]];
      
    // }))/n;

    logger.add_iteration(0,t2.next_time(),msse,0,0,zero_filled_seq,0);

  }



};

#endif //QUANT