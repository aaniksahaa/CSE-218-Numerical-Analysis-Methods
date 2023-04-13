import numpy as np

# Read data from text file into numpy array

file = open("data.txt", "r")  # Change here

separator = " "  # Change here

lines = file.read().split('\n')

xs = np.zeros(0)
ys = np.zeros(0)
values = []
for line in lines:
    line = line.rstrip()
    xy = line.split(separator)
    if(len(xy)==2):
        xs = np.append(xs,float(xy[0]))
        ys = np.append(ys,float(xy[1]))

print(xs)
print(ys)