#Graphing multiple initialization functions, msse by k and time by k
import log_object

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from numpy.polynomial.polynomial import polyfit

def compare_runtimes(log_objects_list,num_diff_methods,klist):

  fig, ax1 = plt.subplots()

  ys = []
  b_totals = []
  m_totals = []


  for i in range(num_diff_methods):
    ys.append([log_objects_list[i][j].total_time for j in range(len(log_objects_list[i]))])
    plt.scatter(klist[0:len(log_objects_list[i])],ys[i],label=log_objects_list[i][0].name)
    b_total, m_total = polyfit(klist[0:len(log_objects_list[i])],ys[i],1)
    b_totals.append(b_total)
    m_totals.append(m_total)
    plt.plot(klist[0:len(log_objects_list[i])],b_totals[i]+m_totals[i]*np.array(klist[0:len(log_objects_list[i])]),"-")

  plt.title(f"Initializations time to k")
  plt.xlabel("k")
  plt.ylabel("total time (s)")


  plt.legend(bbox_to_anchor=(1.2,1.2))
  graph_name = input("Please input name of output graph file scaling k: ")

  plt.tight_layout()

  plt.savefig(graph_name)
  plt.cla()
  plt.clf()
    
  return 0

def compare_msses(log_objects_list,num_diff_methods,klist):

  fig, ax1 = plt.subplots()

  ys = []
  b_totals = []
  m_totals = []

  for i in range(num_diff_methods):
      
    ys.append([log_objects_list[i][j].msse for j in range(len(log_objects_list[i]))])
    plt.scatter(klist[0:len(log_objects_list[i])],ys[i],label=log_objects_list[i][0].name)
    plt.plot(klist[0:len(log_objects_list[i])],ys[i],label=log_objects_list[i][0].name)
    b_total, m_total = polyfit(np.log(klist[0:len(log_objects_list[i])]),ys[i],1)
    b_totals.append(b_total)
    m_totals.append(m_total)
    plt.plot(klist[0:len(log_objects_list[i])],b_totals[i]+m_totals[i]*np.log(np.array(klist[0:len(log_objects_list[i])])),"-")

  plt.title(f"Initialization msse comparison with scaling k")
  plt.xlabel("k")
  plt.ylabel("msse")

  plt.legend(bbox_to_anchor=(1.2,1.2))
  plt.tight_layout()

  graph_name = input("Please input name of output graph file scaling k: ")
  plt.savefig(graph_name)
  plt.cla()
  plt.clf()
    
  return 0

#####Call these graph methods



num_diff_methods = int(input("How many methods? "))
num_ks = int(input("How many ks?"))
klist = []
for i in range(num_ks):
  klist.append(int(input("Input a value of k: ")))


#2D list
log_objects_list = []
for i in range(num_diff_methods):
  log_objects_list.append([])
  num_ks_recorded = int(input("How many ks were recorded for this iter method? "))
  for j in range(num_ks_recorded):
    log_objects_list[i].append(log_object.LogInitObject(input("Please give a csv file name: ").strip()))


compare_runtimes(log_objects_list,num_diff_methods,klist)

compare_msses(log_objects_list,num_diff_methods,klist)