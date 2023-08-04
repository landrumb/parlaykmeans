#file to write random code to pull out interesting facts about the data

import log_object

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from numpy.polynomial.polynomial import polyfit

#input("Please give a csv file name: ").strip()
my_log = log_object.LogObject("/Users/andrewbrady/Documents/GitHub/parlaykmeans/experiments/experiment4/safe_output/test_lazy_naive_100000.csv")

assign_times = [float(my_log.iter_infos[i][1]) for i in range(my_log.niter)]
assign_avg = sum(assign_times)/len(assign_times)
print ("assign avg: ", assign_avg)

