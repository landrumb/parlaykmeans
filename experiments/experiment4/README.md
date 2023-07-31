This is the description of experiment4.

on fern.

Run on this commit:

23d980412f3a45703d27ff08512f9ee6378aaf04

The point of this experiment is to compare the performance of naive and yy, on a bigann slice of 1'000'000 points. \

All of the outputs are placed in this experiment3/output folder, safeguarding files in other folders.
But they should be moved to a safe_output directory to prevent future runs of this experiment from overwriting previous results.

To run this experiment, type
sh experiment3/full_base_experiment.sh 
from parlaykmeans's main directory.

Then you will need to type in Y to confirm the experiment run.

***Caution: in yy there are a couple iterations 'added' in the setup stage (in which logger.add_iter is called)
-- yy still runs 20 true iterations, but the graphs should be shifted leftward by ~3 (should fix in future experiments)
-- the point of this was to show different stages of setup time of yy by itself