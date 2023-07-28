This is the description of experiment7 on fern.

Run on this commit: 0e7aea7b710c6ce125c45fc5629e6ba173d392db

The point of this experiment is to compare the performance of quant to the other methods

All of the outputs are placed in this experiment7/output folder, safeguarding files in other folders.
But they should be moved to a safe_output directory to prevent future runs of this experiment from overwriting previous results.

To run this experiment, type
sh experiment6/full_base_experiment.sh 
from parlaykmeans's main directory.

Then you will need to type in Y to confirm the experiment run.

use

echo "Y" | nohup sh experiment6/full_base_experiment.sh &

__________________________


Remark: I noticed that all of the time was being used in the final centers/assign step.
To confirm whether this time is necessary I added another log add iteration then reran the code 
to get the output in safe_output2.