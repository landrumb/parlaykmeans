This is a test on the LSH initialization on 1billion scale k-means

laxman2@aware:~/parlaykmeans$ sh init_aware.sh
INFO: Analyzed target //:init_test_run (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //:init_test_run up-to-date:
  bazel-bin/init_test_run
INFO: Elapsed time: 15.821s, Critical Path: 15.64s
INFO: 3 processes: 1 internal, 2 local.
INFO: Build completed successfully, 3 total actions
Using Euclidean distance
Detected 1000000000 points with dimension 100
starting bench
Lazy: time-- 11.9759, msse-- 2.6784e+06
n: 1000000000, d: 100, k: 100000, Distance : euclidean
Using Euclidean distance
Detected 1000000000 points with dimension 100
starting bench
starting lsh
Parlay time: Generated the hps: 0.0175
Parlay time: got the hashes: 268.6927
boutaa sort
Parlay time: Just sorted: 5.2836
calculating centers
Finished with center added, now onto asg
Parlay time: Calculated the centers: 96.1766
Parlay time: Got the asg too: 1.7035
left lsh code 
LSH INIT: time-- 373.5, msse-- -nan
n: 1000000000, d: 100, k: 100000, Distance : euclidean

Here's the printout. Was run on init_aware.sh

This experiment shows:
(1) Euclidean distance function buggy randomly?
(2) get_hash, calculate centers steps take up significant amounts of time ... how to improve?
(3) (assuming) LSH does a much better job than Lazy?? (we don't necessarily know thi) 
