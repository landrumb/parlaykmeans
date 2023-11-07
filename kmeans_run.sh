#!bin/bash

#kmeans_test_run is the target for running our main kmeans.cpp program
bazel build kmeans_test_run

#-k 200: set the number of centers to 200
#-i filename: choose this dataset 
#-f bin: related to filetype (always keep set to bin). This is vestigal, previous versions of the code had an option to parse fvecs files, but this version no longer has that feature.
#-t uint8 point data type (choose between uint8, int8, float)
#-D fast: use the EuclideanDistanceFast object for the distance (this is one of the distance objects that uses a vectorized distance function).
#-m 5: max # of iter is 5
#-epsilon 0.01: set epsilon to 0.01 (each iteration we compare epsilon to the max center move to see if we have converged and can stop early)
#-two stable: route the program to the bench_two_stable function, which will do a Lazy initialization (set the centers to be the first k points) then run Naive and Yinyang on that initialization.
./bazel-bin/kmeans_test_run -k 200 -i /ssd1/anndata/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D fast -m 5 -epsilon 0.01 -two stable 
