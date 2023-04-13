import numpy as np

def showMatrix(A):
    print()
    A = np.array(A)
    rows, cols = A.shape
    for i in range(rows):
        for j in range(cols):
            print(f"{A[i][j].round(4):<7}",end="    ")
        print()
    print()

a = [[-1.25463,-23.4567043],[2.356,1.200]]

showMatrix(a)
