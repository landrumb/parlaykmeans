#Brings data files from many csv's and puts them into one csv
#This one for exact k-means methods

import csv
import struct
import log_object

num_diff_methods = int(input("How many methods are we comparing? "))

if (num_diff_methods < 1):
  print("Too few methods, exiting")
  exit(50)


num_ks = int(input("how many diff vals of k are we looking at? "))
klist = []
for i in range(num_ks):
  klist.append(int(input("Input a value of k: ")))

#2D list
log_objects_list = []
for i in range(num_diff_methods):
  log_objects_list.append([])
  for j in range(num_ks):
    #LogObject or LogInitObject depending on use
    log_objects_list[i].append(log_object.LogObject(input("Please give a csv file name: ")))


blocker = input("blocker? ")
if(blocker != "block"):
  exit(3)
output_file_name = input("What csv file shall we output to? ")

csv_data = []
row_data = ["k"]
for rep in range(2):
  for i in range(num_diff_methods):
    row_data.append(log_objects_list[i][0].runner_name)

csv_data.append([item for item in row_data]) #copy

for i in range(num_ks):
  row_data = []
  current_k = klist[i]
  row_data.append(current_k)
  for j in range(num_diff_methods):
    row_data.append(log_objects_list[j][i].total_time)
  for j in range(num_diff_methods):
    row_data.append(log_objects_list[j][i].msse)
  csv_data.append([item for item in row_data]) #copy


with open(output_file_name,'w',newline='\n') as csvfile:
  writer = csv.writer(csvfile)
  writer.writerows(csv_data)