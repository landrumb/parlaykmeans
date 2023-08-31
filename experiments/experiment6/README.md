This is the description of experiment6 on fern.

Run on this commit: 9878ad6f0a7ab9f9960aca48d959a9ebba2e9dd3 
(I think I forgot to commit right before the run so not entirely sure)

The point of this experiment is to compare the performance of naive and yy, on 1 threads to all of the threads available. 

All of the outputs are placed in this experiment6/output folder, safeguarding files in other folders.
But they should be moved to a safe_output directory to prevent future runs of this experiment from overwriting previous results.

To run this experiment, type
sh experiment6/full_base_experiment.sh 
from parlaykmeans's main directory.

Then you will need to type in Y to confirm the experiment run.

use

echo "Y" | nohup sh experiment6/full_base_experiment.sh &


____

EDIT AUGUST 31,2023: 
I added an excel file containing nice graphs of the comparsion between naive and yinyang for scaling k, and 
comparing naive and yy on aware vs fern and for 1 thread vs. 192 threads. The excel file is called 
'Naive yy data graphs.xlsx"