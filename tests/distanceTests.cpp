#include <gtest/gtest.h>
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

TEST(EuclideanDistanceSmall,OneDistance) {
  Distance* D = new EuclideanDistanceSmall();
  int N = 200;
  parlay::sequence<float> ones(N,1);
  parlay::sequence<float> twos(N,2);
  //distance between ones and twos should be the number of elements we are checking, which is i
  for (int i = 0; i <= N; i++) {
    EXPECT_EQ(D->distance(ones.begin(),twos.begin(),i),i);
  }


}

TEST(EuclideanDistance,OnesDistance) {
  Distance* D = new EuclideanDistance();
  int N = 200;
  parlay::sequence<float> ones(N,1);
  parlay::sequence<float> twos(N,2);
  //distance between ones and twos should be the number of elements we are checking, which is i
  for (int i = 8; i <= N; i++) {
    EXPECT_EQ(D->distance(ones.begin(),twos.begin(),i),i);
  }


}