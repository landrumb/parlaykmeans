/* 
    A 'lazy' iteration implementation which does no iterations
 */

#include <tuple>

template <typename T>
struct Lazy {
    void cluster(T* v, size_t n, size_t d, size_t k, T* centers, size_t* asg, size_t max_iter, double epsilon) {
        return;
    }

    void operator()(T* v, size_t n, size_t d, size_t k, float* centers, size_t* asg) {
        return;
    }
};

