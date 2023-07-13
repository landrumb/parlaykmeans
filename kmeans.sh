#!bin/bash
cd ~/parlaykmeans
bazel build kmeans_test_run

P=/ssd1/data
G=/home/landrum/outputs
O=/ssd1/results

BP=$P/bigann
# ./bazel-bin/kmeans_test_run -k 1000 -i $BP/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -m 10

TP=$P/text2image1B
./bazel-bin/kmeans_test_run -k 1000 -i $TP/base.1B.fbin.crop_nb_1000000 -f bin -t float -D Euclidean -m 10

#other implementations
./bazel-bin/kmeans_test_run -k 1000 -i $TP/base.1B.fbin.crop_nb_1000000 -f bin -t float -D Euclidean -m 10 -two yes 
