#include "include/utils/parse_files.h"
//#include "include/utils/parse_command_line.h"
#include "include/utils/NSGDist.h"
#include "lazy.h"
#include "initialization.h"
#include "naive.h"
#include "yinyang_faithful.h"

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"


//T : the data type of a coordinate of our points
//initializer: function used to initialize the centers
//runner: function used to do kmeans (repeated Lloyd's iterations)
//v: flat array of point coordinates
//k: number of clusters requested
//d: dimension of each point
//c: flat array of size k*d that we will write the centers into
//asg: flat array of size n that we will write the point assignments into
//cost: the SSE after the kmeans algorithm finishes
//dist_choice: string representing which distance function we want to use
//max_iter: maximum number of iterations we will run
//epsilon: threshold for comparing the new centers to the old centers; if the 
//  sum of the distance between the new and old centers is less than epsilon, 
//  we stop early
//returns: a double, the cost of the final center choice    
// template<typename T, typename Initializer, typename Runner>         
// void Kmeans(T* v, size_t n, size_t d, size_t k, float* centers, size_t* asg, 
// std::string dist_choice, size_t max_iter = 1000, double epsilon=0.01) {

//     Distance D;
//     if (dist_choice=="euclidean") {
//         if (k >= 36 && d >= 36) {
//             D = EuclideanDistance();

//         }
//         else {
//             D = EuclideanDistanceSmall();
//         }
//     }
//     else {
//         std::cout << "Invalid distance choice" << std::endl;
//         abort();
//     }

//     Kmeans<T,Initializer,Runner>(v,n,d,k,centers,asg,D,max_iter,epsilon);
// }

template<typename T, typename Initializer, typename Runner>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
Distance& D, size_t max_iter = 1000, double epsilon=0.01) {

    Initializer init;
    init(v,n,d,k,c,asg,D);
    Runner run;
    run.cluster(v,n,d,k,c,asg,D,max_iter,epsilon);

}

//debugging function
void debug_dist(Distance& D) {
    parlay::sequence<float> buf1(50,1);
    parlay::sequence<float> buf2(50,2);
    std::cout << "dist init: " << D.distance(make_slice(buf1).begin(),make_slice(buf2).begin(),50) << std::endl;
    std::cout << "finished with dist\n";

    abort(); //cutting early

}


//debugging framework until we set up a proper parsing system
int main() {

    std::cout << "Starting program" << std::endl;

    std::cout << "Note: k, n, d artificially low rn" << std::endl;
    size_t k = 40;
    
    auto file_parts = parse_uint8bin("Data/base.1B.u8bin.crop_nb_1000");
    uint8_t* v = (uint8_t*) std::get<0>(file_parts);
    size_t n = std::get<1>(file_parts);
    size_t d = std::get<2>(file_parts);

    float* c = new float[k*d]; // centers
    size_t* asg = new size_t[n];

    std::cout << "Printing out initial points: " << std::endl;
    for (size_t i = 0; i < n ;i++) {
        for (size_t j = 0; j < d; j++) {
            std::cout << static_cast<int>(v[i*d+j]) << " ";
        }
        std::cout << std::endl;
    }


    std::string dist_choice = "euclidean";
    //DISTANCE MUST BE A dynamically allocated pointer*
    Distance* D;

    if (dist_choice=="euclidean") {
        if (k >= 36 && d >= 36) {
            std::cout << "using vec dist" << std::endl;
            D = new EuclideanDistance();

        }
        else {
            std::cout << "using small dist" <<std::endl;
            D = new EuclideanDistanceSmall();
        }
    }
    else {
        std::cout << "Invalid distance choice" << std::endl;
        abort();
    }

   //debug_dist(*D);

   float* c2 = new float[k*d];
   size_t* asg2 = new size_t[n];

    LazyStart<uint8_t> init;
    init(v,n,d,k,c,asg,*D);
    for (size_t i = 0; i < k*d; i++) {
        c2[i] = c[i];
    }
    for (size_t i = 0; i < n; i++) {
        asg2[i] = asg[i];
    }
    size_t max_iter = 10;
    double epsilon = 0.01;

    NaiveKmeans<uint8_t> nie;
    nie.cluster(v,n,d,k,c,asg,*D,max_iter,epsilon);
    YinyangFaithful<uint8_t> yy;
    yy.cluster(v,n,d,k,c2,asg2,*D,max_iter,epsilon);


    //not actually running kmeans right now
    //Kmeans<uint8_t,LazyStart<uint8_t>,Lazy<uint8_t>>(v,n,d,k,c,asg,D,10,0.01);
    //Kmeans<uint8_t,LazyStart<uint8_t>,NaiveKmeans<uint8_t>>(v,n,d,k,c,asg,*D,10,0.01);
    

    std::cout << "Printing out final centers: "  << std::endl;
    for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < d; j++) {
            std::cout << c[i*d + j] <<  " " << c2[i*d+j] << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    std::cout << "Printing out final assignments: " << std::endl;
    for (size_t i = 0; i < n; i++) {
        std::cout << asg[i] << " " << asg2[i] << std::endl;
    }
    std::cout << std::endl << std::endl;

    delete[] c;
    delete[] asg;
    delete D;


}