/* 
    Initialization methods for k-means clustering.

    methods should be expressed as structs with an operator()

    It's expected that each implementation will:
        - update the centers array with the initial centers
        - update the asg array with the initial assignments

    Implementations which do not update asg (for algorithms which will have to recompute it regardless in the first iteration) should be named with a trailing '_noasg' (e.g. Forgy_noasg)
*/

#include "parlay/random.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "utils/NSGDist.h"
#include "utils/accumulator.h"

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
        // actually uses the lighter weight i % k version
        parlay::parallel_for(0, n, [&](size_t i) {
            asg[i] = i % k;
        });

        //compute centers
        accumulator<float>* acc = new accumulator<float>[k*d];
        parlay::parallel_for(0, n, [&](size_t i) {
            size_t c = i % k;
            for (size_t j = 0; j < d; j++) {
                acc[c*d + j].add(v[i*d + j]);
            }
        });

        parlay::parallel_for(0, k, [&](size_t i) {
            for (size_t j = 0; j < d; j++) {
                size_t count = n / k + (i < n % k);
                centers[i*d + j] = acc[i*d + j].total() / count;
            }
        });
        
        delete[] acc;
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