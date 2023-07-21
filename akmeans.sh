#!bin/bash
bazel build kmeans_test_run

# ./bazel-bin/kmeans_test_run -k 50 -i Data/diff_europe.bin -f bin -t float -D short -m 10 -two yes 


./bazel-bin/kmeans_test_run -k 10 -i Data/europediff_csv.csv -f bin -t float -D short -m 10 -two yes 
