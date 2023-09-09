#Take in a file name, output graphs looking at changes in various variables over iter#
#Small file because meant to be run by a bash script

import log_object

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from numpy.polynomial.polynomial import polyfit


filename = input("Insert csv file name: ")
my_log_obj = log_object.LogObject(filename)
my_log_obj.iter_times_iter_nums_graphs(input("Iter graph file name? "))
my_log_obj.center_reasg_distance_calcs(input("Msse graph file name? "))

