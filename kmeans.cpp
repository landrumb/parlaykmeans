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
#include "quantized.h"

#define INITIALIZER MacQueen
#define INITIALIZER_NAME "MacQueen"
#define RUNNER Naive
#define RUNNER_NAME "Naive"


const float epsilon = 0.00001; //want to stop


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

//bench two kmeans methods on the same data
template <typename T, typename Initializer, typename Runner1, typename Runner2>
inline void bench_two(T* v, size_t n, size_t d, size_t k, Distance& D, 
size_t max_iter=1000, double epsilon=0) {
    std::cout << "Running bench two " << std::endl;

    float* c = new float[k*d]; // centers
    size_t* asg = new size_t[n];

   float* c2 = new float[k*d];
   size_t* asg2 = new size_t[n];
   float* c3 = new float[k*d];
   size_t* asg3 = new size_t[n];

    LazyStart<T> init;
    init(v,n,d,k,c,asg,D);
    for (size_t i = 0; i < k*d; i++) {
        c2[i] = c[i];
        c3[i]=c[i];
    }
    for (size_t i = 0; i < n; i++) {
        asg2[i] = asg[i];
        asg3[i]=asg[i];
    }

    QuantizedKmeans<T> quant;
    kmeans_bench logger_quant = kmeans_bench(n,d,k,max_iter,epsilon,"Lazy","Quant");
    logger_quant.start_time();
    quant.cluster(v,n,d,k,c,asg,D,logger_quant,max_iter,epsilon);
    logger_quant.end_time();

 
    // NaiveKmeans<T> nie;
    // kmeans_bench logger_nie = kmeans_bench(n,d,k,max_iter,
    // epsilon,"Lazy","NaiveKmeans");
    // logger_nie.start_time();
    // nie.cluster(v,n,d,k,c,asg,D,logger_nie,max_iter,epsilon);
    // logger_nie.end_time();
    // std::cout << "cutting out after my naive" << std::endl;
    // abort();
    // std::cout << "starting naive" << std::endl;

    YinyangSimp<T> yy;
    kmeans_bench logger_yy = kmeans_bench(n,d,k,max_iter,epsilon,
    "Lazy","YY");
    logger_yy.start_time();

    yy.cluster(v,n,d,k,c2,asg2,D,logger_yy, max_iter,epsilon);

    logger_yy.end_time();
    std::cout << "Cutting out after yy" << std::endl;
    abort();

    Naive<T> ben_naive;
    kmeans_bench logger = 
    kmeans_bench(n, d, k, max_iter, epsilon, "Lazy", "Naive");
    logger.start_time();
    ben_naive.cluster(v,n,d,k,c3,asg3,D,logger,max_iter,epsilon);
    logger.end_time();

    std::cout << "Printing out first 10 final centers, the first 10 dim: "  << std::endl;
    for (size_t i = 0; i < std::min((size_t) 10,k); i++) {
        for (size_t j = 0; j < std::min((size_t) 10,d); j++) {
            std::cout << c[i*d + j] <<  "|" << c2[i*d+j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Printing out 100 final assignments: " << std::endl;
    for (size_t i = 0; i < std::min(n,(size_t) 100); i++) {
        std::cout << asg[i] << " " << asg2[i] << " " << asg3[i] << std::endl;
        
    }
    std::cout << std::endl << std::endl;

    delete[] c;
    delete[] asg;

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
    std::string use_bench_two = std::string(P.getOptionValue("-two","no"));

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
    } else if (dist=="short") {
        std::cout << "Using short Euclidean" << std::endl;
        D = new EuclideanDistanceSmall();
    } 
    else {
        std::cout << "Error: distance type not specified correctly, specify Euclidean or mips" << std::endl;
        abort();
    }

    if (ft == "bin"){
        if (tp == "float") {
            auto [v, n, d] = parse_fbin(input.c_str());
            if (use_bench_two == "yes") {
                bench_two<float,LazyStart<float>,NaiveKmeans<float>,YinyangSimp<float>>(v,n,d,k,*D,max_iterations,epsilon);

            }
            else {
                bench<float, INITIALIZER<float>, RUNNER<float>>(v, n, d, k, *D, max_iterations, epsilon);


            }
        } else if (tp == "uint8") {
            auto [v, n, d] = parse_uint8bin(input.c_str());
            if (use_bench_two=="yes") {
                bench_two<uint8_t,LazyStart<uint8_t>,NaiveKmeans<uint8_t>,YinyangSimp<uint8_t>>(v,n,d,k,*D,max_iterations,epsilon);
            }
            else {
                bench<uint8_t, INITIALIZER<uint8_t>, RUNNER<uint8_t>>(v, n, d, k, *D, max_iterations, epsilon);


            }
        } else if (tp == "int8") {
            auto [v, n, d] = parse_int8bin(input.c_str());
            if (use_bench_two == "yes") {
                bench_two<int8_t,LazyStart<int8_t>,NaiveKmeans<int8_t>,YinyangSimp<int8_t>>(v,n,d,k,*D,max_iterations,epsilon);

            }
            else {
                bench<int8_t, INITIALIZER<int8_t>, RUNNER<int8_t>>(v, n, d, k, *D, max_iterations, epsilon);


            }
        } else {
            //  this should actually be unreachable
            std::cout << "Error: bin type can only be float, uint8, or int8. Supplied type is " << tp << "." << std::endl;
            abort();
        }
    }

    return 0;

}