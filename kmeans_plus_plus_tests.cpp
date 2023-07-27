//Kmeans++ harder to test but will try

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
#include "include/nisk_kmeans.h"

TEST(KmeansPlusPlus,ValidPointsCheck) {

  auto [v2,n2,d2] = parse_uint8bin("/Users/andrewbrady/Documents/GitHub/parlaykmeans/Data/base.1B.u8bin.crop_nb_1000");//parse_fbin("/ssd1/data/text2image1B/base.1B.fbin.crop_nb_1000");

  //avoid the local reference to binding capture error
  uint8_t* v = v2;
  size_t n = n2;
  size_t d = d2;
  
  size_t k = 100;
  float* c = new float[k*d];
  size_t* asg = new size_t[n];
  Distance* D = new EuclideanDistanceSmall();


  KmeansPlusPlus<uint8_t> init;
  init(v,n,d,k,c,asg,*D);

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t j) {
    return j;
  });

  //confirm all assignments are valid
  for (size_t i = 0; i < n; i++) {
    EXPECT_LE(asg[i],k-1);
    EXPECT_GE(asg[i],0);
  }
  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn, [&] (size_t j) {
    return std::make_pair(asg[j],j);
  }));

  EXPECT_EQ(pts_grouped_by_center.size(),k);

  parlay::sequence<size_t> rang = parlay::tabulate(k,[&] (size_t j) {
    return j;
  });

  size_t* min_asg = new size_t[n];

  //closest point check: as Kmeans++ does assign points to closest center
  parlay::parallel_for(0,n,[&] (size_t i) {
    float buf[2048];
    for (size_t j = 0; j < d; j++) {
      buf[j] = v[i*d+j];
    }
    auto distances = parlay::map(rang, [&] (size_t j) {

      return D->distance(buf,c+j*d,d);
    });
    min_asg[i] = parlay::min_element(distances) - distances.begin();


  });

  for (size_t i = 0; i < n; i++) {
    EXPECT_EQ(min_asg[i],asg[i]);
  }


  delete[] c;
  delete[] asg;
  delete[] min_asg;




}


TEST(Lazy,ValidPointsCheck) {

  auto [v2,n2,d2] = parse_uint8bin("/Users/andrewbrady/Documents/GitHub/parlaykmeans/Data/base.1B.u8bin.crop_nb_1000");//parse_fbin("/ssd1/data/text2image1B/base.1B.fbin.crop_nb_1000");

  //avoid the local reference to binding capture error
  uint8_t* v = v2;
  size_t n = n2;
  size_t d = d2;
  
  size_t k = 100;
  float* c = new float[k*d];
  size_t* asg = new size_t[n];
  Distance* D = new EuclideanDistanceSmall();


  LazyStart<uint8_t> init;
  init(v,n,d,k,c,asg,*D);

  parlay::sequence<size_t> rangn = parlay::tabulate(n,[&] (size_t j) {
    return j;
  });

  //confirm all assignments are valid
  for (size_t i = 0; i < n; i++) {
    EXPECT_LE(asg[i],k-1);
    EXPECT_GE(asg[i],0);
  }
  auto pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn, [&] (size_t j) {
    return std::make_pair(asg[j],j);
  }));

  EXPECT_EQ(pts_grouped_by_center.size(),k);

  delete[] c;
  delete[] asg;


}