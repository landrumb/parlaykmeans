//product quantized k-means (in the style of Jegou's product quantization for nns)
//**RANDOM SEG FAULT BEHAVIOR !!! 
#ifndef PQ_Kmeans
#define PQ_Kmeans

template <typename T>
struct PQKmeans {

  //let's set m=8, k = 2^64
  void cluster(T* v, size_t n, size_t d, size_t k, 
float* c, size_t* asg, Distance& D, kmeans_bench& logger, size_t max_iter, double epsilon,bool suppress_logging=false) {

  if (!suppress_logging) std::cout << "running pq" << std::endl;
  
  if (d > 2048) {
    std::cout << "d greater than 2048, too big, printing d: " << 
    d << std::endl;
    abort();
  }

  parlay::internal::timer t = parlay::internal::timer();
  t.start(); //start timer
  
  if(n == 0) return; //default case 

  size_t m = 8; //note that changing from 8 may cause errors
  size_t kstar = 256;

  if (d % m != 0) {
    std::cout << "need d to be divisible by m, aborting " << std::endl;
    abort();
  }
  size_t smalld = d/m;

  //at index i*kstar*smalld + j*smalld + k, you will find the kth coordinate of the jth center of the center in subvector set i
  float* block_clusters = new float[m*kstar*smalld];
  size_t* block_asg = new size_t[n]; 

  

  //next we quantize each point
  //for each point we want its membership in each of the m=8 blocks
  //at index m*p + j we have the pth point's block id in subvector j
  uint8_t* quantized_points = new uint8_t[m*n];

  //no point in parallelizing because the internal cluster already so parallel? TODO
  for (size_t i = 0; i < m; i++) {
    //copy asg into block_asg; we don't care about the block_asg output
    //but we don't want to overwite asg so we need a copy
    parlay::parallel_for(0,n,[&] (size_t j) {
      block_asg[j]=asg[j];
    
    });

    //we must first initialize the block clusters!
    //we do a Lazy init here
    parlay::parallel_for(0,kstar,[&] (size_t r) {
      for (size_t j = 0; j < smalld; j++) {
        //v[r*d+smalld*i+j] gives the smalld*i+j th coordinate of the rth point in v (smalld*i is our offset for subvector i)
        block_clusters[i*kstar*smalld + r*smalld + j] = v[r*d+smalld*i+j];
      }

    });

    //print init clusters
    // std::cout << "print init clusters " << std::endl;
    // for (size_t r = 0; r < kstar; r++) {
      
    //   for (size_t j = 0; j < smalld; j++) {
    //     std::cout << block_clusters[i*kstar*smalld + r*smalld+j] << " ";
    //   }
    //   std::cout << std::endl;
    // }
    

    kmeans_bench logger_int = kmeans_bench(n,smalld,kstar,5,0,"Outer","Internal PQ Kmeans run");
    //logging suppress set to true TODO, currently set false for debugging
    logger_int.start_time();
    cluster_offset(v,n,smalld,kstar,block_clusters + i*kstar*smalld,block_asg,D,logger_int,5,0,d,smalld*i,false);
    logger_int.end_time();

    //record the blocks each point is assigned to
    parlay::parallel_for(0,n,[&] (size_t j) {
      quantized_points[j*m+i] = block_asg[j];

    });
  }

  std::cout << "quantization took " << t.next_time() << std::endl << std::endl << std::endl;

  //now we proceed to doing (almost normal) kmeans


  size_t iterations = 0;
  float max_diff = 0;
  auto rangk = parlay::iota(k);
  auto rangn = parlay::iota(n);
  float* center_calc_float = new float[k*d]; //do calculations for compute center inside here

  //for each center c, for each block j, we must know how far the subvector c_j of that center is from each of the centers (0 ... kstar-1) in block j
  //at the k*i + m*j + lth index we see the ith center's jth subvector's distance (squared) from the lth block center
  float* center_dists = new float[k*m*kstar];
  float setup_time,assignment_time,update_time=0;

  while (iterations < max_iter) {
    iterations++;

    //we start by precomputing distances from centers to blocks
    //computing in this order as to try to avoid cache misses
    parlay::parallel_for(0,k,[&] (size_t i) {
      for (size_t j = 0; j < m; j++) {
        for (size_t l = 0; l < kstar; l++) {
          center_dists[k*i+m*j+l] = D.distance(c+i*d+smalld*j,block_clusters+kstar*smalld*j+l*smalld,smalld);
        }
      }

    });

    setup_time += t.next_time();

    // Assign each point to the closest center
    parlay::parallel_for(0, n, [&](size_t p) {
    
      
      auto distances = parlay::delayed::map(rangk, [&](size_t i) {
          float my_sum = 0;
          for (size_t j = 0; j < m; j++) {
            my_sum += center_dists[k*i+m*j+quantized_points[p*m+j]];
          }
          return my_sum;
      });

      asg[p] = min_element(distances) - distances.begin();
    });

    assignment_time = t.next_time();

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

    update_time += t.next_time();

    float msse = parlay::reduce(parlay::map(rangn,[&] (size_t i) { 
      float buf[2048];
      T* it = v+i*d;
      for (size_t i = 0; i < d; i++) buf[i]=*(it++);
      return D.distance(buf,c+asg[i]*d,d);
    }))/n; //calculate msse
    
    setup_time += t.next_time(); //setup_time counts msse calculation time
    if (!suppress_logging) logger.add_iteration(assignment_time,update_time,
    msse, 0, 0, deltas,setup_time);
    setup_time=0;
    assignment_time=0;
    update_time=0;

    if (max_diff <= epsilon) break;
    
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
  
  //std::cout << "sup0" << std::endl;
  while (iterations < max_iter) {
    iterations++;

    // std::cout << "Printing centers: " <<std::endl;
    // for (size_t i = 0; i < k; i++) {
    //   for (size_t j = 0; j < d; j++) {
    //     std::cout << c[i*d+j] << " ";
    //   }
    //   std::cout << std::endl;
    // }

    // Assign each point to the closest center
    parlay::parallel_for(0, n, [&](size_t p) {
      float buf[2048];
      T* it = v+p*dorig + offset;
      //note the i-offset!! (because our buf should still start at 0!)
      for (size_t i = 0; i < d; i++) buf[i]=*(it++);
      
      auto distances = parlay::delayed::map(rangk, [&](size_t j) {
          return D.distance(buf, c+j*d,d);
      });

      asg[p] = min_element(distances) - distances.begin();
      //this didn't get caught!
      if (asg[p] > k) {
        std::cout << "bad assign " << p << std::endl;
        for (size_t i = 0; i < distances.size(); i++) {
          std::cout << distances[i] << " ";

        }
        std::cout<<std::endl;
        abort();
      }
    });

    float assignment_time = t.next_time();

      //std::cout << "sup1" << std::endl;


    // Compute new centers
     //copy center coords into center_calc_float
    parlay::parallel_for(0,k*d,[&] (size_t i) {
        center_calc_float[i] = 0;
    });
    //group points by center
    parlay::sequence<std::pair<size_t,parlay::sequence<size_t>>> pts_grouped_by_center = parlay::group_by_key(parlay::map(rangn,[&] (size_t i) {
    return std::pair(asg[i],i);
    }));

    //checking asg:
    for (size_t i = 0; i < n; i++) {
      if (asg[i] > k) {
        std::cout << "bad asg: " << i << " " << asg[i] <<std::endl;
        abort();
      }
    }

     // std::cout << "sup2" << std::endl;

    //add points
    //parlay::parallel_for(0,k,[&] (size_t i) {
    for (size_t i = 0; i < k; i++) { //sequential for debugging
        //  std::cout << "i : " << i << std::endl;
        //  std::cout << "pts_grouped_by_center[i].first: " << pts_grouped_by_center[i].first << std::endl;
        //  std::cout << "list: second size" << pts_grouped_by_center[i].second.size() << std::endl;
        size_t picked_center_d = pts_grouped_by_center[i].first*d;
        for (size_t j = 0; j < pts_grouped_by_center[i].second.size(); j++) {
          size_t point_coord = pts_grouped_by_center[i].second[j]*dorig;
          for (size_t coord = 0; coord < d; coord++) {
            center_calc_float[picked_center_d + coord] += static_cast<float>(v[point_coord + offset+ coord]);
          }
        }
    }
   //},1);
     // std::cout << "sup3" << std::endl;


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
   
    // std::cout << "sup4" << std::endl;

    parlay::sequence<float> deltas = parlay::tabulate(k, [&] (size_t i) {
      return D.distance(center_calc_float+i*d, c + i*d,d);
    });

    max_diff = *parlay::max_element(deltas);
   
    //copy back over centers
    parlay::parallel_for(0,k*d,[&](size_t i) {
      c[i] = center_calc_float[i];
    });

    float update_time = t.next_time();

    float msse = parlay::reduce(parlay::map(rangn,[&] (size_t p) { 
      float buf[2048];
      T* it = v+p*dorig + offset;
      for (size_t i = 0; i < d; i++) buf[i]=*(it++);
      return D.distance(buf,c+asg[p]*d,d);
    }))/n; //calculate msse
    
    float setup_time = t.next_time(); //setup_time counts msse calculation time
    if (!suppress_logging) logger.add_iteration(assignment_time,update_time,
    msse, 0, 0, deltas,setup_time);
      //std::cout << "sup5" << std::endl;


    if (max_diff <= epsilon) break;
    
  }
}

};

#endif //PQKmeans