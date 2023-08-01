#!bin/bash

#file for benching initializers 

bazel build init_test_run

# ./bazel-bin/init_test_run -k 10 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy

# ./bazel-bin/init_test_run -k 10 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH

# ./bazel-bin/init_test_run -k 10000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Lazy

# ./bazel-bin/init_test_run -k 10000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c LSH

# ./bazel-bin/init_test_run -k 10 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Lazy

# ./bazel-bin/init_test_run -k 10 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c LSH

# ./bazel-bin/init_test_run -k 10 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Lazy

# ./bazel-bin/init_test_run -k 10000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c LSH


# ./bazel-bin/init_test_run -k 100000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Lazy

# ./bazel-bin/init_test_run -k 100000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c LSH

./bazel-bin/init_test_run -k 100000 -i /ssd1/data/MSSPACEV1B/spacev1b_base.i8bin -f bin -t uint8 -D short -c Lazy

./bazel-bin/init_test_run -k 100000 -i /ssd1/data/MSSPACEV1B/spacev1b_base.i8bin -f bin -t uint8 -D short -c LSH
