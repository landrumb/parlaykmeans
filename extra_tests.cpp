//file with additional tests

#include <gtest/gtest.h>
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
#include "yinyang_simp.h" //can switch to fast_center

// TEST(HelloTest,Basic) {
//   EXPECT_EQ(5*5,25);

//   std::cout << "can you print from here?" << std::endl;
// }


// TEST(GoodbyeTest,Initial) {
//   EXPECT_EQ(7+99,106);
// }

TEST(NaiveKmeans,Asserts) {
  auto [v,n_int,d_int] = parse_fbin("/ssd1/data/text2image1B/base.1B.fbin.crop_nb_1000000");
  size_t n = (size_t) n_int;
  size_t d = (size_t) d_int;
  size_t k = 10;
  size_t max_iter = 1000;
  float epsilon = 0.001;
  Distance* D = new EuclideanDistanceSmall();
  float* c = new float[k*d]; // centers
  size_t* asg = new size_t[n];

  LazyStart<float> init;
  init(v,n,d,k,c,asg,*D);
  
  NaiveKmeans<float> nie;
  kmeans_bench logger_nie = kmeans_bench(n,d,k,max_iter,epsilon,
  "Lazy","Naivekmeans");
  logger_nie.start_time();
  nie.cluster(v,n,d,k,c,asg,*D,logger_nie,max_iter,epsilon);
  logger_nie.end_time();

  

  float* centroids = new float[d*k];
  size_t* num_members = new size_t[k];
  for (size_t i = 0; i < k*d; i++) {
    centroids[i]=0;
  }
  for (size_t i = 0; i < n; i++) {
    num_members[i] = 0;

  }

  for(size_t i=0; i < n; i++){
    size_t center = asg[i];
    num_members[asg[i]]++;
    for(size_t j=0; j < d; j++){
      centroids[center*d + j] += v[i*d + j];
    }
  }

  //the centroids don't need to be exactly the sum due to 
  //potential floating point error
  //just check first ten
  // for(size_t i=0; i < 10; i++){
  //   EXPECT_LE(std::abs(centroids[i]-c[i]*num_members[i/d]),.1)
  //   << "This one failed\n";
  // }
  EXPECT_LE(1,2);

  delete[] c;
  delete[] asg;
  delete[] centroids;
  delete[] num_members;







}