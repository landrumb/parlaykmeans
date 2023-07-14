//file with additional tests for yy

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

//Test fixture for NaiveKmeans Testing
//Runs clustering on a data file
//The SetUpTestSuite function allows different tests to use the same
//data run
class YyTest : public ::testing::Test {
  protected:

  
  static void SetUpTestSuite()  {
    auto [v2,n_int,d_int] = parse_fbin("/ssd1/data/text2image1B/base.1B.fbin.crop_nb_1000000");
  n = (size_t) n_int;
   d = (size_t) d_int;
   v=v2;
  k = 10;
  max_iter = 1000;
   epsilon = 0; //run to completion
   D = new EuclideanDistanceSmall();

  c = new float[k*d]; // centers
  asg = new size_t[n];

  std::cout << "made it here" << std::endl;

  LazyStart<float> init;
  init(v,n,d,k,c,asg,*D);

  std::cout << "madre it there" << std::endl;
  
  YinyangSimp<float> nie;
  kmeans_bench logger_nie = kmeans_bench(n,d,k,max_iter,epsilon,
  "Lazy","Naivekmeans");
  logger_nie.start_time();
  nie.cluster(v,n,d,k,c,asg,*D,logger_nie,max_iter,epsilon);

  std::cout << "HOLA" << std::endl;

  logger_nie.end_time();


  }

  static void TearDownTestSuite()  {
    delete[] c;
  delete[] asg;

  }
  static float* c;
  static float* v;
  static size_t n;
  static size_t k;
  static size_t d;
  static size_t max_iter;
  static Distance* D;
  static float epsilon;
  static size_t* asg;

};
//initialize static vals
float* YyTest::c = nullptr;
float* YyTest::v = nullptr;
size_t YyTest::n = 0;
size_t YyTest::k = 0;
size_t YyTest::d = 0;
size_t YyTest::max_iter = 0;
Distance* YyTest::D = nullptr;
float YyTest::epsilon = 0;
size_t* YyTest::asg = nullptr;

//Tests that the values in num_members make sense,
//And that points are assigned to a valid center
TEST_F(YyTest,MemberCheck) {

    parlay::sequence<size_t> num_members = parlay::sequence<size_t>(k,0);

   for (size_t i = 0; i < n; i++) {
    num_members[asg[i]]++;

  }

  for (size_t i = 0; i < k; i++) {
    std::cout << "num_mem[" << i << "]: " << num_members[i] << std::endl; 
  }
  std::cout << "total mems " << parlay::reduce(num_members) << std::endl;

  EXPECT_EQ(parlay::reduce(num_members),n);


  //confirm that assignments are to a possible center
  for (size_t i = 0; i < n; i++) {
    EXPECT_LE(asg[i],k-1);
    EXPECT_GE(asg[i],0);
  }

  //confirm that num_members elts are between 0 andn n
  for (size_t i = 0; i < k; i++) {
    EXPECT_GE(num_members[i],0);
    EXPECT_LE(num_members[i],n);
  }

}

//Test on Naive to confirm correct outputs 
//Tests whether the center is indeed the center of the points
TEST_F(YyTest,CenterCentroids) {
  
  double* centroids = new double[d*k];
  parlay::sequence<size_t> num_members = parlay::sequence<size_t>(k,0);
  for (size_t i = 0; i < k*d; i++) {
    centroids[i]=0;
  }

  parlay::parallel_for(0,d,[&] (size_t j) {
    for (size_t i = 0; i < n; i++) {
      size_t center = asg[i];
      
      centroids[center*d+j] += v[i*d+j];
    }

  });

  for (size_t i = 0; i < n; i++) {
    num_members[asg[i]]++;

  }

  for (size_t i = 0; i < k; i++) {
    std::cout << "num_mem[" << i << "]: " << num_members[i] << std::endl; 
  }
  std::cout << "total mems " << parlay::reduce(num_members) << std::endl;

  EXPECT_EQ(parlay::reduce(num_members),n);

  //confirm that assignments are to a possible center
  for (size_t i = 0; i < n; i++) {
    EXPECT_LE(asg[i],k-1);
    EXPECT_GE(asg[i],0);
  }

  //confirm that num_members elts are between 0 andn n
  for (size_t i = 0; i < k; i++) {
    EXPECT_GE(num_members[i],0);
    EXPECT_LE(num_members[i],n);
  }

  std::cout << "HERE4" << std::endl;

  for (size_t i = 0; i < 10; i++) {
    std::cout << centroids[i] << " " << c[i] << std::endl;
  }

  //TODO how tightly should we set this?
  //putting centroids in double to get a bit better precision
  for(size_t i=0; i < k*d; i++){
    EXPECT_LE(std::abs(centroids[i]/num_members[i/d]-c[i]),.1)
    << "This one failed\n";
  }
  
  delete[] centroids;

 
   
}

//Makes sure that each point is assigned to its closest center
TEST_F(YyTest,ClosestPoints) {
   //rang is range from 0 to k incl excl
    parlay::sequence<size_t> rang = parlay::tabulate(k,[&] (size_t i) {
      return i;
    });

    size_t* bests = new size_t[n];

    parlay::parallel_for(0,n,[&] (size_t i) {
      auto distances = parlay::delayed::map(rang, [&](size_t j) {
            return D->distance(v+i*d, c+j*d,d); });

     bests[i] = min_element(distances) - distances.begin();
    });

    for (size_t i = 0; i < n; i++) {
      EXPECT_EQ(bests[i],asg[i]);
    }

    delete[] bests; //memory cleanup
}