import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

x = []
y = []
z = []
xc = []
yc = []
zc = []
cluster = []
n = 0
k = 0
d = 0
# read file

#colors = ['r', 'g', 'b', 'y', 'm', 'c', 'k', ]


# Open the CSV file
with open('data.csv', 'r') as file:
    # Create a CSV reader object
    csv_reader = csv.reader(file)

    itr = 0
    # Read each row of the CSV file
    for row in csv_reader:
        if itr == 0:
            n = int(row[0])
            k = int(row[1])
            d = int(row[2])
            colors = random.sample(range(k), k)
        else:    
            if itr < n+1: 
                # Access the data in each field of the row
                for i in range(0,d+1):
                    if i == 0:  
                        x.append(float(row[i]))
                    elif i == 1 and d > 1:
                        y.append(float(row[i]))
                    elif i == 2 and d > 2:
                        z.append(float(row[i]))
                    elif i == d:
                        cluster.append(colors[int(row[i])])
            else:
                for i in range(0,d):
                     if i == 0:  
                        xc.append(float(row[i]))
                     elif i == 1 and d > 1:
                        yc.append(float(row[i]))
                     elif i == 2 and d > 2:
                        zc.append(float(row[i]))


        itr = itr + 1
            
#print(x)
#print(y)
#print(z)
#print(cluster)

print(d)

if d > 2:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c=cluster, marker='o')
    ax.scatter(xc, yc, zc, c=colors, marker='x')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Scatter Plot')
elif d > 1: 
    plt.scatter(x, y, c=cluster, marker ='o')
    plt.scatter(xc, yc, c= colors, marker ='x')
    # Set plot title and labels
    plt.title('2D Scatter Plot with Colors')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')



# Display the plot
plt.savefig('scatter_plot.png')
plt.show()