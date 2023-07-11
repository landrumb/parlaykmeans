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

#define INITIALIZER MacQueen
#define INITIALIZER_NAME "MacQueen"
#define RUNNER Naive
#define RUNNER_NAME "Naive"


const float epsilon = 0.0;


template<typename T, typename Initializer, typename Runner>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
Distance& D, kmeans_bench logger, size_t max_iter = 1000, double epsilon=0.01) {

    Initializer init;
    init(v,n,d,k,c,asg,D);

    std::cout << "Initialization complete" << std::endl;

    Runner run;
    run.cluster(v,n,d,k,c,asg,D,logger,max_iter,epsilon);

    std::cout << "Clustering complete" << std::endl;

}

template <typename T, typename initializer, typename runner>
inline void bench(T* v, size_t n, size_t d, size_t k, Distance& D, size_t max_iter = 1000, double epsilon=0) {
    float* centers = new float[k*d];
    size_t* asg = new size_t[n];
    kmeans_bench logger = kmeans_bench(n, d, k, max_iter, epsilon, INITIALIZER_NAME, RUNNER_NAME);
    logger.start_time();
    Kmeans<T, initializer, runner>(v, n, d, k, centers, asg, D, logger, max_iter, epsilon);
    logger.end_time();

    return;
}

int main(int argc, char* argv[]){
    commandLine P(argc, argv, "[-k <n_clusters>] [-m <iterations>] [-o <output>] [-i <input>] [-f <ft>] [-t <tp>] [-D <dist>]");

    size_t k = P.getOptionLongValue("-k", 10); // k is number of clusters
    size_t max_iterations = P.getOptionLongValue("-m", 1000); // max_iterations is the max # of Lloyd iters kmeans will run
    std::string output = std::string(P.getOptionValue("-o", "kmeans_results.csv")); // maybe the kmeans results get written into this csv
    std::string input = std::string(P.getOptionValue("-i", "")); // the input file
    std::string ft = std::string(P.getOptionValue("-f", "bin")); // file type, bin or vecs
    std::string tp = std::string(P.getOptionValue("-t", "uint8")); // data type
    std::string dist = std::string(P.getOptionValue("-D", "Euclidian")); // distance choice

    if(input == ""){ // if no input file given, quit
        std::cout << "Error: input file not specified" << std::endl;
        abort();
    }

    if((ft != "bin") && (ft != "vec")){ // if the file type chosen is not one of the two approved file types, quit 
    std::cout << "Error: file type not specified correctly, specify bin or vec" << std::endl;
    abort();
    }

    if((tp != "uint8") && (tp != "int8") && (tp != "float")){ // if the data type isn't one of the three approved data types, quit
        std::cout << "Error: vector type not specified correctly, specify int8, uint8, or float" << std::endl;
        abort();
    }

    if((ft == "vec") && (tp == "int8")){ // you can't store int8s in a vec file apparently I guess
        std::cout << "Error: incompatible file and vector types" << std::endl;
        abort();
    }

    // TODO: add support for vec files
    if (ft == "vec") {
        std::cout << "Error: vec file type not supported yet" << std::endl;
        abort();
    }

    Distance* D; // create a distance object, it can either by Euclidian or MIPS
    if (dist == "Euclidean") { 
        std::cout << "Using Euclidean distance" << std::endl;
        D = new EuclideanDistance();
    } else if (dist == "mips") {
        std::cout << "Using MIPS distance" << std::endl;
        D = new Mips_Distance();
    } else {
        std::cout << "Error: distance type not specified correctly, specify Euclidean or mips" << std::endl;
        abort();
    }

    if (ft == "bin"){
        if (tp == "float") {
            auto [v, n, d] = parse_fbin(input.c_str());
            bench<float, INITIALIZER<float>, RUNNER<float>>(v, n, d, k, *D, max_iterations, epsilon);
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            bench<uint8_t, INITIALIZER<uint8_t>, RUNNER<uint8_t>>(v, n, d, k, *D, max_iterations, epsilon);
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            bench<int8_t, INITIALIZER<int8_t>, RUNNER<int8_t>>(v, n, d, k, *D, max_iterations, epsilon);
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}