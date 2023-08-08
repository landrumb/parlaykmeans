At long last, an experiment to measure the performance of bench centers.
Update centers consists of two steps (in one view): 
a group by to see which points belong to which centers,
and
an add, where the points are added to the center they belong to.

Our group by is the standard parlay group_by_key currently.

We run two different methods of adding -
one that adds with a reduce parallelised over k and d
one that adds in sequence, just parallelised over k, to reduce
cache misses ideally

The purpose is to see if we currently have a sufficiently fast method for update_centers on the scale of interest, n=1bil, k=100K.

Run on fern starting August 8.

Please note: this experiment was run before changing the granularity on add_points2 to 1 (my bad should have committed before changing this): so if this experiment is rerun the times for add2 will be even faster than the current data for k=10.

Run on this commit:
ea2b791d6e7e9387494cad284a02a731b7cebf79

To get *Similar* results (not the exact same due to above issue)