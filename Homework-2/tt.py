import numpy as np

a = np.array([[1,2,3],[4,5,6],[7,8,9]])

b = np.array([[6],[7],[8]])

arr = np.hstack((a,b))

print(arr)

print(arr.shape[1])

for i in range(5,-1,-1):
    print(i)

a[1] = a[1] - 3
print(a)

c = np.array(a)

print(c)