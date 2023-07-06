/* 
    Initialization methods for k-means clustering.

    methods should be expressed as structs with an operator()
*/

#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "NSGDist.h"
#include "accumulator.h"

/* 
    Randomly assigns each point to a cluster, then computes the centers of each 
    cluster.
 */
template<typename T>
struct Forgy {
    /* 
    args:
        v: pointer to flat array of points
        n: number of points
        k: number of clusters
        d: dimension of points
        centers: pointer to array of centers
        asg: pointer to array of assignments
     */
    void operator()(T* v, size_t n, size_t d, size_t k, float* centers, 
size_t* asg) {
        parlay::random r;
        parlay::parallel_for(0, n, [&](size_t i) {
            asg[i] = r.ith_rand(i) % k;
        });

        parlay::group
    }
};


//Lazy start makes the first k points the first k centers
//Then assigns cyclically
template<typename T>
struct LazyStart {
    void initialize(T* v, size_t n, size_t d, size_t k, float* c, 
    size_t* asg, Distance& D) {
        for (size_t i = 0; i < k; i++) {
            for (size_t j = 0; j < d; j++) {
                c[i*d + j] = v[i*d + j];

            }

        }
        for (size_t i = 0; i < n; i++) {
            asg[i] = i % k;
        }
    }
};