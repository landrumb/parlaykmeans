//product quantized k-means (in the style of Jegou's product quantization for nns)
#ifndef PQ_Kmeans
#define PQ_Kmeans

template <typename T>
struct PQKmeans {

  //let's set m=8, k = 2^64
  void cluster(T* v, size_t n, size_t d, size_t k, 
float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon,bool suppress_logging=false) {
  size_t m = 8; //note that changing from 8 may cause errors
  size_t kstar = 256;
  //at index i*kstar*d + j*d + k, you will find the kth coordinate of the jth center of the center in subvector set i
  float* block_clusters = new float[m*kstar*d];

  if (d % m != 0) {
    std::cout << "need d to be divisible by m, aborting " << std::endl;
    abort();
  }
  size_t smalld = d/m;
  for (size_t i = 0; i < m; i++) {
    //continue here//
  }



}

//cluster the given points, starting at a certain offset (as to grab a portion of the coordinates)
//use dorig to get the proper place in the data point v
void cluster_offset(T* v, size_t n, size_t d, size_t k, 
float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon,size_t dorig, size_t offset, bool suppress_logging=false) {

  if (!suppress_logging) std::cout << "running pq" << std::endl;
  
  if (d > 2048) {
    std::cout << "d greater than 2048, too big, printing d: " << 
    d << std::endl;
    abort();
  }

  parlay::internal::timer t = parlay::internal::timer();
  t.start(); //start timer
  
  if(n == 0) return; //default case 

  size_t iterations = 0;
  float max_diff = 0;
  auto rangk = parlay::iota(k);
  auto rangn = parlay::iota(n);
  float* center_calc_float = new float[k*d]; //do calculations for compute center inside here

  while (iterations < max_iter) {
    iterations++;

    // Assign each point to the closest center
    parlay::parallel_for(0, n, [&](size_t i) {
      float buf[2048];
      T* it = v+i*dorig;
      for (size_t i = offset; i < d; i++) buf[i]=*(it++);
      
      auto distances = parlay::delayed::map(rangk, [&](size_t i) {
          return D.distance(buf, c+i*d,d);
      });

      asg[i] = min_element(distances) - distances.begin();
    });

    float assignment_time = t.next_time();

    // Compute new centers
     //copy center coords into center_calc_float
    parlay::parallel_for(0,k*d,[&] (size_t i) {
        center_calc_float[i] = 0;
    });
    //group points by center
    parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
    }));
    //add points
    parlay::parallel_for(0,k,[&] (size_t i) {
        size_t picked_center_d = pts_grouped_by_center[i].first*d;
        for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
          size_t point_coord = pts_grouped_by_center[i].second[j]*d;
          for (size_t coord = 0; coord < d; coord++) {
            center_calc_float[picked_center_d + coord] += static_cast<float>(v[point_coord + coord]);
          }
        }
    },1);

    parlay::parallel_for(0,k,[&] (size_t i) {

      parlay::parallel_for(0,d,[&] (size_t coord) {
        if (pts_grouped_by_center[i].second.size() > 0) {
          center_calc_float[pts_grouped_by_center[i].first*d+coord] /= pts_grouped_by_center[i].second.size();
        }
        else { //if no points belong to this center
          center_calc_float[pts_grouped_by_center[i].first*d+coord] = c[pts_grouped_by_center[i].first*d+coord];
        }
      });
    
    });
   
    parlay::sequence<float> deltas = parlay::tabulate(k, [&] (size_t i) {
      return D.distance(center_calc_float+i*d, c + i*d,d);
    });

    max_diff = *parlay::max_element(deltas);
   
    //copy back over centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      c[i] = center_calc_float[i];
    });

    float update_time = t.next_time();

    float msse = parlay::reduce(parlay::map(rangn,[&] (size_t i) { 
      float buf[2048];
      T* it = v+i*dorig;
      for (size_t i = offset; i < d; i++) buf[i]=*(it++);
      return D.distance(buf,c+asg[i]*d,d);
    }))/n; //calculate msse
    
    float setup_time = t.next_time(); //setup_time counts msse calculation time
    if (!suppress_logging) logger.add_iteration(assignment_time,update_time,
    msse, 0, 0, deltas,setup_time);

    if (max_diff <= epsilon) break;
    
  }
}

};

#endif //PQKmeans