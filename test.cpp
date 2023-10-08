// #include <stdio.h>
// #include <fstream>
// #include <ostream>
// #include <iostream>

//modifying file for testing git branching

// #include <gtest/gtest.h>

// #include "include/yinyang_simp.h"

// #include "include/utils/parse_files.h"
// #include "include/utils/parse_command_line.h"
// #include "include/utils/NSGDist.h"
// #include "lazy.h"
// #include "initialization.h"
// #include "naive.h"
// #include "yinyang_simp.h"
// #include "include/utils/kmeans_bench.h"

// #include "parlay/sequence.h"
// #include "parlay/parallel.h"
// #include "parlay/primitives.h"


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

// #include <cmath> 


int main() {
  std::cout << "hello" << std::endl;
  std::cout << parlay::num_workers() << std::endl;

  // sequence<size_t> test;
  // test.push_back(5);
  // test.push_back(22);
  // auto dist = parlay::map(test,[&](size_t j))
  // std::cout << "min of " << std::min(5,7) << std::endl;
  // size_t a = 0;
  // size_t b = 1;
  // std::cout << "min of " << std::min(a,b) << std::endl;

  // std::cout << "Working on parlay filter " << std::endl;

  // parlay::sequence<size_t> my_ints = {5,4,20,10};
  // auto filtered = parlay::filter(my_ints,[&] (size_t& i) {
  //   return i > 7;

  // });
  // //std::cout << filtered << std::endl;
  // for (size_t i = 0; i < filtered.size(); i++) {
  //   std::cout << filtered[i] << " ";
  // }
  //   std::cout << std::endl;


  // parlay::sequence<size_t> mapped = parlay::map(filtered,
  // [&] (size_t i) {
  //   return i+1;
  // });

  // for (size_t i = 0; i < mapped.size(); i++) {
  //   std::cout << mapped[i] << " ";
  // }
  // std::cout << std::endl;

  // size_t k = 81;
  // size_t kstar = 3;
  // size_t split = 4;
  

  // std::cout << "Debugging center selects" << std::endl;

  // std::cout << "pow(Kstar, 0) " << std::pow(kstar,0) << std::endl;
  // std::cout << "1 / pow(kstar,0)" << (1/std::pow(kstar,0)) << std::endl;

  //  for (size_t i = 0; i < k; i++) {
  //     parlay::sequence<size_t> center_selects(split);
  //     std::cout << "printing center select for i: " << i << std::endl;

  //     for (size_t j = 0; j < split; j++) {
  //       center_selects[j] = static_cast<long>(i / std::pow(kstar,j)) % kstar;
  //       std::cout << center_selects[j] << " ";

  //     }
  //     std::cout << std::endl;
    
  //   }
  // parlay::parallel_for(0,20,[&] (size_t i) {
  //   if (i == 10) {
  //     return;

  //   }
  //   std::cout << i << std::endl;

  // });

  std::cout << "using delayed seq " << std::endl;
  auto s1 = parlay::delayed_tabulate(10,[&] (size_t i) {
    return i;
  });
  parlay::delayed_sequence<int, int, int > s2 = parlay::iota(20);
  // for (int i = 0; i < s2.size(); i++) {
  //   std::cout << s2[i] << " ";
  // }
  std::cout << std::endl;

  return 1;
}
