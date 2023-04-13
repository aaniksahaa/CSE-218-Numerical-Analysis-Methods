import numpy as np
import math
import matplotlib.pyplot as plt
import csv

# n is the order of interpolation
# xs should be sorted in ascending order
def lagrangeinterpolate(xs,ys,n,givenx):  

    xs = np.array(xs)
    ys = np.array(ys)

    # Sort the data points according to x values 

    tuples = []
    for i in range(xs.shape[0]):
        tuples.append((xs[i],ys[i]))
    tuples.sort()
    for i in range(xs.shape[0]):
        xs[i] = tuples[i][0]
        ys[i] = tuples[i][1]

    N = xs.shape[0]

    # Checking if enough data points are available
    if(N < (n+1)):
        print("Sorry, Not enough data points for {}th order interpolation".format(n))
        return -1

    # Check if sizes of xs and ys are equal
    if(N != ys.shape[0]):
        print("Sorry, Invalid Data. Size of arrays not equal.")
        return -1

    # Checking if the given x is in between given data points
    if(givenx > xs.max() or givenx < xs.min()):
        print("Sorry, Point to be interpolated out of data range.")
        return -1

    newxs = np.zeros(n+1)
    newys = np.zeros(n+1)

    cnt = 0
    leftx = -1
    rightx = -1

    for i in range(N):
        if(xs[i] == givenx):
            print("Already Given. Interpolated y = {}".format(ys[i]))
        if(xs[i]>givenx):
            newxs[cnt] = xs[i-1]
            newys[cnt] = ys[i-1]
            cnt += 1
            newxs[cnt] = xs[i]
            newys[cnt] = ys[i]
            cnt += 1
            break
    
    # Choosing (n+1) closest data points
    j = i-2
    k = i+1

    while(cnt < n+1):
        if( k>N-1 or (abs(xs[j]-givenx) < abs(xs[k]-givenx))):
            newxs[cnt] = xs[j]
            newys[cnt] = ys[j]
            j-=1
            cnt+=1
        else:
            newxs[cnt] = xs[k]
            newys[cnt] = ys[k]
            k+=1
            cnt+=1

    #print(newxs)
    #print(newys)

    allxs = xs
    allys = ys

    xs = newxs  # Final (n+1) data points
    ys = newys

    result = 0.0
    for i in range(n+1):

        weight = 1.0
        for j in range(n+1):
            if(j != i):
                weight = weight*(givenx - xs[j])
        for j in range(n+1):
            if(j != i):
                weight = weight/(xs[i] - xs[j])
        
        result = result + weight*ys[i]

    return result

def lagrange(xs,ys,maxn,givenx):

    gap = " "
    spaces = 4*" " + "|" + 4*" "
    print("\n\nLagrange Interpolation\n")
    print(100*"-")
    print(5*gap + "Interpolation Order" + spaces + "Given x" + spaces + "Interpolated y" + spaces + "Relative Approx. Error")
    print(100*"-")
    prevy = -1

    for i in range(1,maxn+1):
        interpolated_y = lagrangeinterpolate(xs,ys,i,givenx)
        error = -1
        if(i>1):
            error = (abs((interpolated_y - prevy)/interpolated_y))*100

        firstpart = 14*gap + str(i) + 9*gap + spaces + "%.4f"%givenx + spaces + 3*gap + "%.4f"%interpolated_y + 3*gap + spaces + 6*gap
        
        if(interpolated_y != -1): # if there has been no error
            if(i==1):
                print(firstpart + 8*"-")
            else:
                print(firstpart + "%.4f"%error + " %")

        prevy = interpolated_y
    
    print("\nInterpolated y for {}th order polynomial = {}\n\n".format(maxn,interpolated_y))

N = 100

divdiff = np.zeros((N,N))

# Calculate all divided differences by dynamic programming
def precalc_div_differences(xs,ys):
    xs = np.array(xs)
    ys = np.array(ys)
    n = xs.shape[0]
    for i in range(n):
        for j in range(n-i):
            if(i==0):
                divdiff[j][j+i] = ys[j]
            else:
                divdiff[j][j+i] = (divdiff[j+1][j+i] - divdiff[j][j+i-1])/(xs[j+i]-xs[j])


def newtoninterpolate(xs,ys,n,givenx):   # For a given order

    xs = np.array(xs)
    ys = np.array(ys)

    # Sort the data points according to x values 

    tuples = []
    for i in range(xs.shape[0]):
        tuples.append((xs[i],ys[i]))
    tuples.sort()
    for i in range(xs.shape[0]):
        xs[i] = tuples[i][0]
        ys[i] = tuples[i][1]

    N = xs.shape[0]

    # Checking if enough data points are available
    if(N < (n+1)):
        print("Sorry, Not enough data points for {}th order interpolation".format(n))
        return -1

    # Check if sizes of xs and ys are equal
    if(N != ys.shape[0]):
        print("Sorry, Invalid Data. Size of arrays not equal.")
        return -1

    # Checking if the given x is in between given data points
    if(givenx > xs.max() or givenx < xs.min()):
        print("Sorry, Point to be interpolated out of data range.")
        return -1

    newxs = np.zeros(n+1)
    newys = np.zeros(n+1)

    cnt = 0
    leftx = -1
    rightx = -1

    for i in range(N):
        if(xs[i] == givenx):
            print("Already Given. Interpolated y = {}".format(ys[i]))
        if(xs[i]>givenx):
            newxs[cnt] = xs[i-1]
            newys[cnt] = ys[i-1]
            cnt += 1
            newxs[cnt] = xs[i]
            newys[cnt] = ys[i]
            cnt += 1
            break
    
    # Choosing (n+1) closest data points
    j = i-2
    k = i+1

    while(cnt < n+1):
        if( k>N-1 or (abs(xs[j]-givenx) < abs(xs[k]-givenx))):
            newxs[cnt] = xs[j]
            newys[cnt] = ys[j]
            j-=1
            cnt+=1
        else:
            newxs[cnt] = xs[k]
            newys[cnt] = ys[k]
            k+=1
            cnt+=1

    #print(newxs)
    #print(newys)

    allxs = xs
    allys = ys

    xs = newxs  # Final (n+1) data points
    ys = newys

    # Sort the final data points according to x values 

    tuples = []
    for i in range(xs.shape[0]):
        tuples.append((xs[i],ys[i]))
    tuples.sort()
    for i in range(xs.shape[0]):
        xs[i] = tuples[i][0]
        ys[i] = tuples[i][1]

    #print(xs)
    #print(ys)

    precalc_div_differences(xs,ys)

    result = 0.0
    for i in range(n+1):
        temp = divdiff[0][i]
        for j in range(i):   # from 0 to i-1
            temp = temp*(givenx - xs[j])
        result = result + temp

    return result


def newton(xs,ys,maxn,givenx):   # For all orders upto maxn

    gap = " "
    spaces = 4*" " + "|" + 4*" "
    print("\n\nNewton Interpolation\n")
    print(100*"-")
    print(5*gap + "Interpolation Order" + spaces + "Given x" + spaces + "Interpolated y" + spaces + "Relative Approx. Error")
    print(100*"-")
    prevy = -1

    for i in range(1,maxn+1):
        interpolated_y = newtoninterpolate(xs,ys,i,givenx)
        error = -1
        if(i>1):
            error = (abs((interpolated_y - prevy)/interpolated_y))*100

        firstpart = 14*gap + str(i) + 9*gap + spaces + "%.4f"%givenx + spaces + 3*gap + "%.4f"%interpolated_y + 3*gap + spaces + 6*gap
        
        if(interpolated_y != -1): # if there has been no error
            if(i==1):
                print(firstpart + 8*"-")
            else:
                print(firstpart + "%.4f"%error + " %")

        prevy = interpolated_y
    
    print("\nInterpolated y for {}th order polynomial = {}\n\n".format(maxn,interpolated_y))

# Reading data points from csv file

datax = []
datay = []
dataz = []

def read_data_from_csv(filename):    # reads columns one by one into the lists
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)
    i = 0
    for row in data: 
        i += 1
        if(i==1):     # reads excluding first row containing header
            continue
        datax.append(float(row[0]))
        datay.append(float(row[1]))
        dataz.append(float(row[2]))
        i+=1
    #print(datax)
    #print(datay)

# Drawing Line plot of Given data points and interpolated points 
# Change as per problem       
def linegraph():  
    # sort if needed
    xval = 16    # given x value for which y is unknown
    gx = [xval]
    gy = [newtoninterpolate(datax,datay,4,xval)]
    plt.plot(datax,datay,"o-g",label="P = 1 bar")
    plt.plot(gx,gy,"*-b",label="Interpolated P = 1 bar",markersize=8)
    plt.xlabel('Temparature')
    plt.ylabel('Solutbility')
    #plt.xlim([-0.15,0.2])
    plt.legend(loc="best")
    plt.grid()
    plt.show()

# Test Data

xs = [0,10,15,20,22.5,30]
ys = [0,227.04,362.78,517.35,602.97,901.67]

lagrange(xs,ys,5,16)
newton(xs,ys,5,16)

#read_data_from_csv("dd.csv")
#linegraph()


