#Given different implementations,
#Take in logging data for both,
#Graph the total runtime over scaling k
#Graph the msses over scaling k
#Fix # of iters

import log_object

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from numpy.polynomial.polynomial import polyfit

def compare_runtimes(log_objects_list,num_diff_methods,klist,num_ks):

  fig, ax1 = plt.subplots()

  ys = []

  for i in range(num_diff_methods):
    ys.append([log_objects_list[i][j].total_time for j in range(num_ks)])
    plt.scatter(klist,ys[i],label=log_objects_list[i][0].runner_name)
  plt.title(f"Time comparison with scaling k, num iters {log_objects_list[0][0].max_iter}")
  plt.xlabel("k")
  plt.ylabel("iter time (s)")

  plt.legend()
  graph_name = input("Please input name of output graph file scaling k: ")
  plt.savefig(graph_name)
  plt.cla()
  plt.clf()
    
  return 0

def compare_msses(log_objects_list,num_diff_methods,klist,num_ks):

  fig, ax1 = plt.subplots()

  ys = []

  for i in range(num_diff_methods):
    ys.append([log_objects_list[i][j].msse for j in range(num_ks)])
    plt.scatter(klist,ys[i],label=log_objects_list[i][0].runner_name)
  plt.title(f"Msse comparison with scaling k, num iters {log_objects_list[0][0].max_iter}")
  plt.xlabel("k")
  plt.ylabel("msse")

  plt.legend()
  graph_name = input("Please input name of output graph file scaling k: ")
  plt.savefig(graph_name)
  plt.cla()
  plt.clf()
    
  return 0

num_diff_methods = int(input("How many methods are we comparing? "))

if (num_diff_methods < 1):
  print("Too few methods, exiting")
  exit(50)

init_names = []
run_names = []
for i in range(num_diff_methods):
  init_names.append(input("Please give the name of the initialization method: "))
  run_names.append(input("Please give the name of the runner method: "))



num_ks = int(input("how many diff vals of k are we looking at? "))
klist = []
for i in range(num_ks):
  klist.append(int(input("Input a value of k: ")))

#2D list
log_objects_list = []
for i in range(num_diff_methods):
  log_objects_list.append([])
  for j in range(num_ks):
    log_objects_list[i].append(log_object.LogObject(input("Please give a csv file name: ")))

#Debugging, confirm read in was mostly correct
#Necessary but not sufficient check
for i in range(num_diff_methods):
  for j in range(num_ks):
    if (log_objects_list[i][j].k != klist[j]):
      print("k does not align, exiting")
      exit(61)
    if (log_objects_list[i][j].init_name != init_names[i]):
      print("init name does not align, exiting")
      print(log_objects_list[i][j].init_name, init_names[i])
      exit(62)
    if (log_objects_list[i][j].runner_name != run_names[i]):
      print("run name does not align, exiting")
      exit(63)
  

compare_runtimes(log_objects_list,num_diff_methods,klist,num_ks)

compare_msses(log_objects_list,num_diff_methods,klist,num_ks)