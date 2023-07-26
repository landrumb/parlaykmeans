#Log object for a kmeans run, that holds the necessary info
#in a sensible way

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
from numpy.polynomial.polynomial import polyfit


class LogObject:

  def __init__(self,filename):

    with open (filename,"r") as data_file:
      csv_reader = csv.reader(data_file)

      itr = 0

      #note that next removes the row from the csv reader
      overall_info_names = next(csv_reader)


      overall_info = next(csv_reader)
      self.n = int(overall_info[0])
      self.d = int(overall_info[1])
      self.k = int(overall_info[2])
      self.max_iter = int(overall_info[3])
      self.epsilon = float(overall_info[4])
      self.niter = int(overall_info[5])
      self.msse = float(overall_info[6])
      self.total_time = float(overall_info[7])
      self.init_name = overall_info[8].strip()
      self.runner_name = overall_info[9].strip()

      per_iter_names = next(csv_reader)

      self.iter_infos = []
      for row in csv_reader:
        self.iter_infos.append(row)


  def iter_times_iter_nums_graphs(self,graph_name):
    #first graph: iter time (y) to niter (n)

    iter_times = []
    for row in self.iter_infos:
      iter_times.append(float(row[1]) + float(row[2]) + float(row[3]))

    iter_asg_times = []

    for row in self.iter_infos:
      iter_asg_times.append(float(row[1]))

    iter_update_times = []

    for row in self.iter_infos:
      iter_update_times.append(float(row[2]))

    iter_setup_times = []

    for row in self.iter_infos:
      iter_setup_times.append(float(row[3]))

    iter_nums = []
    for row in self.iter_infos:
      iter_nums.append(int(row[0]))

    print("iter_nums.size()", len(iter_nums))
    print("iter times.size()", len(iter_times))

    fig = plt.figure()
    ax = fig.add_subplot(111)


    plt.scatter(iter_nums,iter_times,label="total times")
    plt.scatter(iter_nums,iter_asg_times,label="asg times")
    plt.scatter(iter_nums,iter_update_times,label="upd times")
    plt.scatter(iter_nums,iter_setup_times,label="setup times")

    #line of best fit for total times
    b_total, m_total = polyfit(iter_nums,iter_times,1)
    plt.plot(iter_nums,b_total+m_total*np.array(iter_nums),"-")

    b_asg, m_asg = polyfit(iter_nums,iter_asg_times,1)
    plt.plot(iter_nums,b_asg+m_asg*np.array(iter_nums),"-")

    b_update, m_update = polyfit(iter_nums,iter_update_times,1)
    plt.plot(iter_nums,b_update+m_update*np.array(iter_nums),"-")

    b_setup, m_setup = polyfit(iter_nums,iter_setup_times,1)
    plt.plot(iter_nums,b_setup+m_setup*np.array(iter_nums),"-")

    plt.title("Iter times to iter#")
    plt.xlabel("iter#")

    plt.ylabel("iter time (s)")


    plt.legend()

    plt.savefig(graph_name)

    #plt.show()

    plt.cla()
    plt.clf()


  #TODO use a dictionary to get which row is which stat, instead of hardcoding the #s
  #this would be useful in case the stats change over time
  def center_reasg_distance_calcs(self,graph_name):
    iter_nums = [int(row[0]) for row in self.iter_infos]
    distance_calcs = [int(row[5]) for row in self.iter_infos]
    center_reasgs = [int(row[6]) for row in self.iter_infos]

    fig, ax1 = plt.subplots()

    color1 = 'tab:red'
    ax1.set_xlabel("iter#")
    ax1.set_ylabel("distance calcs",color=color1)
    p1 = ax1.plot(iter_nums,distance_calcs,color=color1,label = "distance calcs")

    ax2 = ax1.twinx()
    color2 = 'tab:blue'
    ax2.set_ylabel("center reasg", color=color2)
    p2 = ax2.plot(iter_nums,center_reasgs,color=color2,label="center reasgs")
    fig.tight_layout() #helps stuff show

    center_move_means = [float(row[7]) for row in self.iter_infos]
    ax3 = ax1.twinx()
    color3 = 'tab:green'
    ax3.set_ylabel("center move means", color=color3)
    p3 = ax3.plot(iter_nums,center_move_means,color=color3,label="center mean moves")
    ax3.spines['right'].set_position(('outward', 60))

    msses = [float(row[4]) for row in self.iter_infos]
    ax4 = ax1.twinx()
    color4 = 'tab:orange'
    ax4.set_ylabel("msse", color=color4)
    p4 = ax4.plot(iter_nums,msses,color=color4,label="msses")
    ax4.spines['right'].set_position(('outward', 120))

    #plt.legend()
    ax1.legend(handles=p1+p2+p3+p4,loc='best')

    plt.savefig(graph_name,bbox_inches='tight')
    plt.cla()
    plt.clf()
