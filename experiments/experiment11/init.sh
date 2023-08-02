#!bin/bash

#file for benching initializers 

bazel build init_test_run

./bazel-bin/init_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c MacQueen -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Forgy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c KmeansPlusPlus -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 100 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 100 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c MacQueen -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 100 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Forgy -o "./experiments/experiment11/output/"

 ./bazel-bin/init_test_run -k 100 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c KmeansPlusPlus -o "./experiments/experiment11/output/"

 ./bazel-bin/init_test_run -k 100 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH -o "./experiments/experiment11/output/"


./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c MacQueen -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Forgy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH -o "./experiments/experiment11/output/"

# # ./bazel-bin/init_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c KmeansPlusPlus


./bazel-bin/init_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c MacQueen -o "./experiments/experiment11/output/"

./bazel-bin/init_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Forgy -o "./experiments/experiment11/output/"


./bazel-bin/init_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH -o "./experiments/experiment11/output/"

