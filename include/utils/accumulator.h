#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

#include <math.h>
#include <stdlib.h>

// no performance difference between 64 and higher values, but 32 is significantly slower
#define BUFFER 64 

template <typename T>
struct accumulator {

    size_t width = BUFFER / sizeof(T);
    parlay::sequence<T> counts;

    accumulator() {
        counts = parlay::sequence<T>(parlay::num_workers() * width, (T)0);
    }

    void increment(){
        counts[parlay::worker_id()*width]++;
    }

    void decrement(){
        // probably not really valid on unsigned types but as long as the total is positive overflow should in theory make this work as intented anyways
        counts[parlay::worker_id()*width]--;
    }

    void add(T val){
        counts[parlay::worker_id()*width] += val;
    }

    void subtract(T val){
        counts[parlay::worker_id()*width] -= val;
    }

    T total(){
        return parlay::reduce(counts);
    }

    void reset(){
        counts = parlay::sequence<T>(parlay::num_workers() * width, (T)0);
    }
};

template <typename T>
struct minimizer {

    size_t width = BUFFER / sizeof(T);
    parlay::sequence<T> counts;

    minimizer(T init) {
        counts = parlay::sequence<T>(parlay::num_workers() * width, (T)init);
    }

    void update(T val){
        counts[parlay::worker_id()*width] = std::min(counts[parlay::worker_id()*width], val);
    }

    T total(){
        return parlay::reduce(counts, parlay::minm<T>());
    }

    void reset(){
        counts = parlay::sequence<T>(parlay::num_workers() * width, (T)0);
    }
};

template <typename T>
struct maximizer {

        size_t width = BUFFER / sizeof(T);
        parlay::sequence<T> counts;

        maximizer(T init) {
            counts = parlay::sequence<T>(parlay::num_workers() * width, (T)init);
        }
    
        void update(T val){
            counts[parlay::worker_id()*width] = std::max(counts[parlay::worker_id()*width], val);
        }
    
        T total(){
            return parlay::reduce(counts, parlay::maxm<T>());
        }
    
        void reset(){
            counts = parlay::sequence<T>(parlay::num_workers() * width, (T)0);
        }
};
