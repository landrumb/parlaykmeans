This experiment was a simple calculation of the average assign time
for the naive algorithm on n=1mil, k=100K. 

Using the data_scripter.py program, we take the assign time data
from
parlaykmeans/experiments/experiment4/safe_output/test_lazy_naive_100000.csv
and find the average.

This results in a average of 

(base) andrewbrady@guest-wireless-10-106-110-65 parlaykmeans % python includ
e/utils/data_scripter.py
assign avg:  217.77204999999998
(base) andrewbrady@guest-wireless-10-106-110-65 parlaykmeans % 

218s. 

Thus for n=1bil I would expect that we would see an assign runtime of
218'000s, or 2.5 days. 
Note that this 218s result was from a fern run (experiment 4).

To reproduce this assign average, run on commit:
 245adacfe68ae66dbbe2ba4373763bacac77ef98

 

