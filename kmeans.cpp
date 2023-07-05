
//T : the data type of a coordinate of our points
//initializer: function used to initialize the centers
//runner: function used to do kmeans (repeated Lloyd's iterations)
//v: flat array of point coordinates
//k: number of clusters requested
//d: dimension of each point
//centers: flat array of size k*d that we will write the centers into
//asg: flat array of size n that we will write the point assignments into
//cost: the SSE after the kmeans algorithm finishes
//dist_choice: string representing which distance function we want to use
//max_iter: maximum number of iterations we will run
//epsilon: threshold for comparing the new centers to the old centers; if the sum of the distance between the new and old centers is less than epsilon, we stop early
//returns: a double, the cost of the final center choice
template<typename T, typename initializer<T>, typename runner<T>>         
void Kmeans(T* v, size_t n, size_t k, size_t d, float* centers, size_t* asg, double& cost, std::string dist_choice, size_t max_iter = 1000, double epsilon=0.01) {

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

    initializer.assign(v,n,k,d,centers,asg,D);
    runner.cluster(v,n,k,d,centers,asg,cost,D,max_iter,epsilon);
    
}