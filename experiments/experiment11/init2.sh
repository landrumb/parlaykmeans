#!bin/bash

bazel build init_test_run

./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c KmeansPlusPlus -o "./experiments/experiment11/output/"
