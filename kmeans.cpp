
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
template<typename T, typename initializer<T>, typename runner<T>>         
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

    Kmeans(v,n,d,k,centers,asg,D,max_iter,epsilon);

}

template<typename T, typename initializer<T>, typename runner<T>>         
void Kmeans(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
Distance& dist, size_t max_iter = 1000, double epsilon=0.01) {

    initializer.assign(v,n,d,k,c,asg,D);
    runner.cluster(v,n,d,k,c,asg,D,max_iter,epsilon);

}


//debugging framework until we set up a proper parsing system
int main() {

    uint8_t* v;
    size_t n;
    size_t d;
    tie(v,n,d) = parse_uint8bin("../Data/base.1B.u8bin.crop_nb_1000");
    float* c = new float[k*d]; // centers
    size_t* asg = new size_t[n];

    std::cout << "Printing out initial points: " << std::endl;
    for (int i = 0; i < n ;i++) {
        for (int j = 0; j < d; j++) {
            std::cout << v[i*d+j] << " ";
        }
        std::cout << std::endl;
    }

    //Kmeans(v,n,d,k,c,asg,"euclidean",10,0.01);

    std::cout << "Printing out final centers: "  << std::endl;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < d; j++) {
            std::cout << c[i*d + j] <<  " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Printing out final assignments: " << std::endl;
    for (int i = 0; i < 1000; i++) {
        std::cout << asg[i] << " ";
    }
    std::cout << std::endl;



}