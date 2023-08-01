This is the description of experiment5.

on aware.

Run on this commit:
1bfe69f9ac0cbf2e2c5fa1ea6fa56d3e19126b44




The point of this experiment is to compare the performance of aware to fern to naive and yy, on a bigann slice of 1'000'000 points.

Just running k=10K (not k=100K) to get runs to finish in reasonable time.

All of the outputs are placed in this experiment5/output folder, safeguarding files in other folders.
But they should be moved to a safe_output directory to prevent future runs of this experiment from overwriting previous results.

To run this experiment, type
sh experiment5/full_base_experiment.sh 
from parlaykmeans's main directory.

Then you will need to type in Y to confirm the experiment run.

I recommend running with nohup as follows:

echo "Y" | nohup sh experiment5/full_base_experiment.sh & 