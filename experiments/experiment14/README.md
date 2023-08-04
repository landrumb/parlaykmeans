This experiment was meant to confirm that our closest point
check does experimentally scale with n linearly

We ran a closest point check for n=10mil to 100mil with 10mil 
increments, k=100, from the spacev1b_base.i8bin dataset

Run with sh assign.sh

Then for redundancy we do the same thing with the 100mil 
slice of BigANN.

Run with experiments/experiment14/assign2.sh

Run on this commit:


andrew@fern:~/parlaykmeans$ sh assign.sh
INFO: Analyzed target //:assign_bench_run (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //:assign_bench_run up-to-date:
  bazel-bin/assign_bench_run
INFO: Elapsed time: 0.134s, Critical Path: 0.01s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              10000000
d:              100
k:              100
total time:     0.125351
assign time:    0.386667
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              20000000
d:              100
k:              100
total time:     0.240107
assign time:    0.760011
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              30000000
d:              100
k:              100
total time:     0.354987
assign time:    1.14582
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              40000000
d:              100
k:              100
total time:     0.473343
assign time:    1.52409
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              50000000
d:              100
k:              100
total time:     0.586191
assign time:    1.90582
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              60000000
d:              100
k:              100
total time:     0.698451
assign time:    2.29348
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              70000000
d:              100
k:              100
total time:     0.828035
assign time:    2.76456
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              80000000
d:              100
k:              100
total time:     0.929165
assign time:    3.15899
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              90000000
d:              100
k:              100
total time:     1.03977
assign time:    3.56124
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              100000000
d:              100
k:              100
total time:     1.13994
assign time:    3.81005
andrew@fern:~/parlaykmeans$ 
_____________________________________


With running the BigANN slice:

andrew@fern:~/parlaykmeans$ sh assign.sh
INFO: Analyzed target //:assign_bench_run (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //:assign_bench_run up-to-date:
  bazel-bin/assign_bench_run
INFO: Elapsed time: 0.134s, Critical Path: 0.01s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              10000000
d:              100
k:              100
total time:     0.125351
assign time:    0.386667
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              20000000
d:              100
k:              100
total time:     0.240107
assign time:    0.760011
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              30000000
d:              100
k:              100
total time:     0.354987
assign time:    1.14582
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              40000000
d:              100
k:              100
total time:     0.473343
assign time:    1.52409
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              50000000
d:              100
k:              100
total time:     0.586191
assign time:    1.90582
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              60000000
d:              100
k:              100
total time:     0.698451
assign time:    2.29348
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              70000000
d:              100
k:              100
total time:     0.828035
assign time:    2.76456
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              80000000
d:              100
k:              100
total time:     0.929165
assign time:    3.15899
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              90000000
d:              100
k:              100
total time:     1.03977
assign time:    3.56124
Using Euclidean distance
Detected 1000000000 points with dimension 100
 initialization with Lazy
n:              100000000
d:              100
k:              100
total time:     1.13994
assign time:    3.81005
andrew@fern:~/parlaykmeans$ 
andrew@fern:~/parlaykmeans$ sh experiments/experiment14/assign2.sh
INFO: Analyzed target //:assign_bench_run (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //:assign_bench_run up-to-date:
  bazel-bin/assign_bench_run
INFO: Elapsed time: 0.135s, Critical Path: 0.01s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              10000000
d:              128
k:              100
total time:     0.13684
assign time:    0.528362
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              20000000
d:              128
k:              100
total time:     0.241031
assign time:    0.930682
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              30000000
d:              128
k:              100
total time:     0.361131
assign time:    1.32129
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              40000000
d:              128
k:              100
total time:     0.479721
assign time:    1.78096
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              50000000
d:              128
k:              100
total time:     0.585799
assign time:    2.19806
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              60000000
d:              128
k:              100
total time:     0.694685
assign time:    2.65247
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              70000000
d:              128
k:              100
total time:     0.8139
assign time:    3.083
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              80000000
d:              128
k:              100
total time:     0.930696
assign time:    3.50628
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              90000000
d:              128
k:              100
total time:     1.04078
assign time:    3.93146
Using Euclidean distance
Detected 100000000 points with dimension 128
 initialization with Lazy
n:              100000000
d:              128
k:              100
total time:     1.1514
assign time:    4.35144
andrew@fern:~/parlaykmeans$ 