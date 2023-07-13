#!bin/bash
bazel build kmeans_test_run

./bazel-bin/kmeans_test_run -k 50 -i Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 10 -two yes 


./bazel-bin/kmeans_test_run -k 10 -i Data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -m 10 -two yes 
