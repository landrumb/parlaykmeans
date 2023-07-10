/* 
    A 'lazy' iteration implementation which does no iterations
 */

#include <tuple>
#include "utils/NSGDist.h"
#include "utils/kmeans_bench.h"

template <typename T>
struct Lazy {
    void cluster(T* v, size_t n, size_t d, size_t k, T* centers, size_t* asg, size_t max_iter, double epsilon, kmeans_bench& logger) {
        return;
    }

    void cluster(T* v, size_t n, size_t d, size_t k, T* centers, size_t* asg, Distance& D, size_t max_iter, double epsilon, kmeans_bench& logger) {
        return;
    }

    void cluster(T* v, size_t n, size_t d, size_t k, float* centers, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
        return;
    }

    void operator()(T* v, size_t n, size_t d, size_t k, float* centers, size_t* asg, kmeans_bench& logger) {
        return;
    }
};

