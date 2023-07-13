#include <stdio.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include <gtest/gtest.h>

#include "include/yinyang_simp.h"

#include "include/utils/parse_files.h"
#include "include/utils/parse_command_line.h"
#include "include/utils/NSGDist.h"
#include "lazy.h"
#include "initialization.h"
#include "naive.h"
#include "yinyang_simp.h"
#include "include/utils/kmeans_bench.h"

#include "include/parlay/sequence.h"
#include "include/parlay/parallel.h"
#include "include/parlay/primitives.h"


// struct kmeans_test : testing::test{

//   size_t k;
//   uint8_t* v;
//   size_t n;
//   size d;
//   float* c;
//   size_t* asg;
//   LazyStart<uint8_t> init;
//   size_t max_iter;
//   double epsilon;
//   YinyangSimp<uint8_t> yy;

//   kmeans_test(){

//       k = 40;

//       auto file_parts = parse_uint8bin("Data/base.1B.u8bin.crop_nb_1000");
//       v = (uint8_t*) std::get<0>(file_parts);

//       n = std::get<1>(file_parts);
//       d = std::get<2>(file_parts);

//       c = new float[k*d]; // centers
//       asg = new size_t[n];

//       std::string dist_choice = "euclidean";
//       Distance* D;
//       D = new EuclideanDistanceSmall();  

//       init(v,n,d,k,c,asg,*D);

//       max_iter = 10;
//       epsilon = 0.01;

//       yy.cluster(v,n,d,k,c2,asg2,*D,max_iter,epsilon);

//     }
// };

// TEST_F(kmeans_test, correct_centroids) {

//   float centroids[d*k];

//   for(size_t i; i < n; i++){
//     size_t center = asg[i];
//     for(size_t j; j < d; j++){
//       centroids[center*d + j] += v[i*d + j]
//     }
//   }

//   for(size t i; i < k*d; i++){
//     ASSERT_EQ(centroids[i], c[i])
//   }
// }

int main() {
  std::cout << "hello" << std::endl;
  std::cout << parlay::num_workers() << std::endl;

  // sequence<size_t> test;
  // test.push_back(5);
  // test.push_back(22);
  // auto dist = parlay::map(test,[&](size_t j))
  return 1;
}
