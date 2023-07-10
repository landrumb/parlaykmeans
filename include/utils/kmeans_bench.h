#ifndef KMEANS_BENCH_H
#define KMEANS_BENCH_H

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/internal/get_time.h"

#include <iostream>

/* 
    stores benchmarking properties for a single iteration of kmeans
*/
struct iteration_bench {
    double assign_time;
    double update_time;
    double msse;
    size_t distance_calculations;
    size_t center_reassignments;
    parlay::sequence<double> center_movements;

    iteration_bench(double assign_time, double update_time, double msse, size_t distance_calculations, size_t center_reassignments, parlay::sequence<double> center_movements) : assign_time(assign_time), update_time(update_time), msse(msse), distance_calculations(distance_calculations), center_reassignments(center_reassignments) {
        this->center_movements = center_movements;
    };

    void print() {
        std::cout << "assignment time:       \t" << assign_time << std::endl;
        std::cout << "update time:           \t" << update_time << std::endl;
        std::cout << "msse:                  \t" << msse << std::endl;
        std::cout << "distance calculations: \t" << distance_calculations << std::endl;
        std::cout << "center reassignments:  \t" << center_reassignments << std::endl;
        std::cout << "center movements:" << std::endl;
        std::cout << "\tmean: " << parlay::reduce(center_movements) / center_movements.size() << std::endl;
        std::cout << "\tmin:  " << parlay::reduce(center_movements, parlay::minm<double>()) << std::endl;
        std::cout << "\tmax:  " << parlay::reduce(center_movements, parlay::maxm<double>()) << std::endl;
        // could also throw quartiles in here
    }
};

/* 
    stores benchmarking properties for a kmeans run
*/
struct kmeans_bench {
    size_t n;
    size_t d;
    size_t k;
    size_t max_iter;
    double epsilon;
    parlay::sequence<iteration_bench> iterations;
    parlay::internal::timer t;
    parlay::internal::timer iteration_timer;
    double total_time = 0.0;
    std::string initializer;
    std::string runner;

    size_t n_iterations = 0;


    kmeans_bench(size_t n, size_t d, size_t k, size_t max_iter, double epsilon, std::string initializer, std::string runner) : n(n), d(d), k(k), max_iter(max_iter), epsilon(epsilon), initializer(initializer), runner(runner){
        iterations = parlay::sequence<iteration_bench>();
    }

    void start_time() {
        std::cout << initializer << " initialization with " << runner << " iterations." << std::endl;
        std::cout << "n:         \t" << n << std::endl;
        std::cout << "d:         \t" << d << std::endl;
        std::cout << "k:         \t" << k << std::endl;
        std::cout << "max_iter:  \t" << max_iter << std::endl;
        std::cout << "epsilon:   \t" << epsilon << std::endl;

        t.start();
        iteration_timer.start();
    }

    void end_time() {
        total_time = t.stop();

        std::cout << "iterations:\t" << n_iterations << std::endl;
        std::cout << "msse:      \t" << iterations[iterations.size() - 1].msse << std::endl;
        std::cout << "total time:\t" << total_time << std::endl;
        std::cout << "avg iteration time:\t" << total_time / iterations.size() << std::endl;
        // std::cout << "total assignment time:\t" << parlay::reduce(parlay::tabulate(iterations, iterations.size(), )) << std::endl;
        // std::cout << "total update time:\t" << parlay::reduce(iterations, parlay::addm<double>(), [](iteration_bench b) {return b.update_time;}) << std::endl;
    }

    /* 
    args:
        assign_time: time to assign points to centers
        update_time: time to update centers
        msse: mean sum squared error
        distance_calculations: number of distance calculations
        center_reassignments: number of points assigned to a new center
        center_movements: sequence of distances moved by each center in an iteration
     */
    void add_iteration(double assign_time, double update_time, double msse, size_t distance_calculations, size_t center_reassignments, parlay::sequence<double> center_movements) {
        iterations.push_back(iteration_bench(assign_time, update_time, msse, distance_calculations, center_reassignments, center_movements));
        n_iterations++;

        std::cout << "iteration " << n_iterations << " complete. (" << iteration_timer.next_time() << ")" << std::endl;
    }

    void print() {
        std::cout << "n:         \t" << n << std::endl;
        std::cout << "d:         \t" << d << std::endl;
        std::cout << "k:         \t" << k << std::endl;
        std::cout << "max_iter:  \t" << max_iter << std::endl;
        std::cout << "epsilon:   \t" << epsilon << std::endl;
        std::cout << "iterations:\t" << n_iterations << std::endl;
        std::cout << "msse:      \t" << iterations[iterations.size() - 1].msse << std::endl;
        std::cout << "total time:\t" << total_time << std::endl;
        std::cout << "avg iteration time:\t" << total_time / iterations.size() << std::endl;
        // std::cout << "total assignment time:\t" << parlay::reduce(parlay::tabulate(iterations, iterations.size(), )) << std::endl;
        // std::cout << "total update time:\t" << parlay::reduce(iterations, parlay::addm<double>(), [](iteration_bench b) {return b.update_time;}) << std::endl;
    }

};

#endif