#!bin/bash

#file for benching initializers 

bazel build update_bench_run


./bazel-bin/update_bench_run -k 10 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -v normal -n 1000000


# ./bazel-bin/update_bench_run -k 10 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000 -v add_only


# ./bazel-bin/update_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000 -v add_only

# ./bazel-bin/update_bench_run -k 1000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000 -v add_only


# ./bazel-bin/update_bench_run -k 10000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000 -v add_only


# ./bazel-bin/update_bench_run -k 10 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000 -v add_only


# ./bazel-bin/update_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000 -v add_only

# ./bazel-bin/update_bench_run -k 1000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000 -v add_only


# ./bazel-bin/update_bench_run -k 10000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000 -v add_only



# ./bazel-bin/update_bench_run -k 1000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/update_bench_run -k 1000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_10000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/update_bench_run -k 1000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/update_bench_run -k 10000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/update_bench_run -k 10000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_10000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/update_bench_run -k 10000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Standard


# ./bazel-bin/assign_bench_run -k 1000 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -c Standard

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_10000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_10000000 -f bin -t uint8 -D short -c Standard

# ./bazel-bin/assign_bench_run -k 10 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/assign_bench_run -k 10 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D short -c Standard

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/bigann/bigann/base.1B.u8bin.crop_nb_100000000 -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/assign_bench_run -k 10 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard

# ./bazel-bin/assign_bench_run -k 10 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D short -c Standard

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 20000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 30000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 40000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 50000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 60000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 70000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 80000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 90000000

# ./bazel-bin/assign_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000