#!bin/bash

#file for benching initializers 

bazel build group_by_bench_run

./bazel-bin/group_by_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000

./bazel-bin/group_by_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 20000000



./bazel-bin/group_by_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 40000000


./bazel-bin/group_by_bench_run -k 100 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 80000000



./bazel-bin/group_by_bench_run -k 1000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 1000000


./bazel-bin/group_by_bench_run -k 10000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 1000000

./bazel-bin/group_by_bench_run -k 1000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000


./bazel-bin/group_by_bench_run -k 10000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000

./bazel-bin/group_by_bench_run -k 1000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 40000000


./bazel-bin/group_by_bench_run -k 10000 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 40000000

./bazel-bin/group_by_bench_run -k 3 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000


./bazel-bin/group_by_bench_run -k 3 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 10000000


./bazel-bin/group_by_bench_run -k 3 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000


./bazel-bin/group_by_bench_run -k 3 -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -n 100000000