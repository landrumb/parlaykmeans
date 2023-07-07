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

#include "include/initialization.h"
#include "include/lazy.h"
#include "include/naive.h"


template<typename T, typename Initializer, typename Runner>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
Distance& D, kmeans_bench logger, size_t max_iter = 1000, double epsilon=0.01) {

    Initializer init;
    init.initialize(v,n,d,k,c,asg,D);
    Runner run;
    run.cluster(v,n,d,k,c,asg,D,max_iter,epsilon);

}

template <typename T, typename initializer, typename runner>
inline void bench(T* v, size_t n, size_t d, size_t k, Distance& D, size_t max_iter = 1000, double epsilon=0.01) {
    float* centers = new float[k*d];
    size_t* asg = new size_t[n];
    kmeans_bench logger = kmeans_bench(n, d, k, max_iter, epsilon);
    logger.start_time();
    Kmeans<T, initializer, runner>(v, n, d, k, centers, asg, D, logger, max_iter, epsilon);
    logger.end_time();
    logger.print();

    // std::cout << empty_centers(centers) << " empty centers" << std::endl;
    // std::cout << avg_center_size(centers) << " average center size" << std::endl;
    // std::cout << dist_between_first_center_and_centroid(v, centers, *D) << " distance between first center and centroid" << std::endl;
    // std::pair<int, int> min_max = min_max_cluster_size(centers);
    // std::cout << min_max.first << " minimum cluster size, " << min_max.second << " maximum cluster size" << std::endl;
    // std::cout << msse(v, centers, *D) << " mean sum squared error" << std::endl;

    return;
}