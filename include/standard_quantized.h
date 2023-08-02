//the other quantized.h file is doing a cartesian k-means approach
//whereas this version is using PQ (via k-means on smaller blocks) to reduce the number
//of dimensions we are looking at

struct StandardQuantizedKmeans {

  void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon) {
    size_t kstar = 64; //number of centers per split
    size_t split = 4; //number of splits we have
  }
}