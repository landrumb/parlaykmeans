This is the description of experiment16.

on fern.

Run on this commit:



The point of this experiment is to compare the performance of naive and yy, on a bigann slice of 1'000'000 points. \

But unlike experiment 4, run on k=1K,2K,3125,4K,5K

All of the outputs are placed in this experiment3/output folder, safeguarding files in other folders.
But they should be moved to a safe_output directory to prevent future runs of this experiment from overwriting previous results.

Remark: even though the code says "Lazy" an LSH initialization was used, this was a typo in the run (note the lower initial msses as well).

