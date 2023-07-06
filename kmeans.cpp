#include "utils/parse_files.h"
#include "lazy.h"
#include "initialization.h"
#include "NSGDist.h"

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
template<typename T, typename Initializer, typename Runner>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* centers, size_t* asg, 
std::string dist_choice, size_t max_iter = 1000, double epsilon=0.01) {

    Distance D;
    if (dist_choice=="euclidean") {
        if (k >= 36 && d >= 36) {
            D = EuclideanDistance();

        }
        else {
            D = EuclideanDistanceSmall();
        }
    }
    else {
        std::cout << "Invalid distance choice" << std::endl;
        abort();
    }

    Kmeans<T,Initializer,Runner>(v,n,d,k,centers,asg,D,max_iter,epsilon);

}

template<typename T, typename Initializer, typename Runner>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
Distance& D, size_t max_iter = 1000, double epsilon=0.01) {

    Initializer init;
    init.initialize(v,n,d,k,c,asg,D);
    Runner run;
    run.cluster(v,n,d,k,c,asg,D,max_iter,epsilon);

}


//debugging framework until we set up a proper parsing system
int main() {

    std::cout << "Starting program" << std::endl;

    size_t k = 10;
    
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

    //not actually running kmeans right now
    Kmeans<uint8_t,LazyStart<uint8_t>,Lazy<uint8_t>>(v,n,d,k,c,asg,"euclidean",10,0.01);

    std::cout << "Printing out final centers: "  << std::endl;
    for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < d; j++) {
            std::cout << c[i*d + j] <<  " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Printing out final assignments: " << std::endl;
    for (size_t i = 0; i < n; i++) {
        std::cout << asg[i] << " ";
    }
    std::cout << std::endl;

    delete[] c;
    delete[] asg;



}