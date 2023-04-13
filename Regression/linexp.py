import numpy as np
import math
import matplotlib.pyplot as plt

def bisection(lo,hi,err,mx,f):
    eps = 10**(-18) 
    cnt = 0
    e = 100 # dummy
    ans = 0.0
    while(e - err > (-eps)):
        prelo = lo
        prehi = hi
        if(cnt == mx):
            break
        mid = (lo+hi)/2
        if(mid == 0.0):
            mid += 0.00001
        #print( "lo = "+str(lo) + "  hi = "+str(hi) + "  mid = "+str(mid))
        fl = f(lo)
        fm = f(mid)
        check = fl*fm
        if(check - 0.0 < eps):  # handling precision error
            hi = mid
        elif(check - 0.0 > -eps):
            lo = mid
        elif(math.fabs(check - 0.0) < eps):
            break
        cnt = cnt + 1
        if(cnt > 1):
            e = (math.fabs((mid-ans)/mid))*100
        ans = mid

    ans = mid
    return ans

def swapRows(A,i,j):
    A[[i,j]] = A[[j,i]]

def showMatrix(A):
    print()
    A = np.array(A)
    rows, cols = A.shape
    for i in range(rows):
        print("      ",end="")
        for j in range(cols):
            print(f"{A[i][j].round(4):<7}",end="    ")
        print()
    print()

def  GaussianElimination(A, B, pivot = 1, showall = 1 ):

    # Convert matrices to numpy array

    A = np.array(A)
    B = np.array(B)

    # Setting printing precision

    np.set_printoptions(precision=4)

    totalswaps = 0

    n = A.shape[0]

    # Formation of Augmented Matrix

    C = np.hstack((A,B))  

    if(showall==1):
        print("Augmented Matrix: ")
        showMatrix(C)

    for i in range(n-1):

        if(showall == 1):
            print("Step {}:\n".format(i+1))

        goodrow = i
        currmax = abs(C[goodrow][i])

        if(pivot == 1):

            # Finding the row with maximum absolute value of intended column
            # This is called partial pivoting

            for j in range(i+1,n):
                if(abs(C[j][i]) > currmax):
                    goodrow = j
                    currmax = C[j][i]
            
            if(currmax == 0):
                print("No nonzero value found")
                print("Sorry, the given equations cannot be solved...\n")
                S = np.zeros((1,2))  # Will be diffentiated by shape
                return S

            # Swapping rows

            if(i != goodrow):
                if(showall == 1):
                    print("   Rows {} and {} are swapped for partial pivoting.".format(i+1,goodrow+1))
                swapRows(C,i,goodrow)
                totalswaps += 1
                if(showall == 1):
                    showMatrix(C)
            
            # Making desired elements zero

            for j in range(i+1,n):
                multiplier = C[j][i]/C[i][i]
                C[j] = C[j] - multiplier*C[i]
                if(showall == 1):
                    print("   {} multiplied by row {} is subtracted from row {}".format(round(multiplier,4),i+1,j+1))
                    showMatrix(C)
        
        else:

            # Avoiding Zero-division error

            if(C[i][i] == 0):

                if(showall == 1):
                    print("   First element is zero. Finding any non-zero element...")
                
                for j in range(i+1,n):
                    if(C[j][i] != 0):
                        goodrow = j
                        currmax = C[j][i]
                        print("   Non-zero found at Row {}\n".format(j+1))
                        if(showall == 1):
                            print("   Rows {} and {} are swapped to avoid zero-division error".format(i+1,goodrow+1))
                        swapRows(C,i,goodrow)
                        totalswaps += 1
                        if(showall == 1):
                            showMatrix(C)
                        break
                
                if(C[goodrow][i]==0):
                    print("No nonzero value found")
                    print("Sorry, the given equations cannot be solved...\n")
                    S = np.zeros((1,2))  # Will be diffentiated by shape
                    return S

            # Making desired elements zero

            for j in range(i+1,n):
                multiplier = C[j][i]/C[i][i]
                C[j] = C[j] - multiplier*C[i]
                if(showall == 1):
                    print("   {} multiplied by row {} is subtracted from row {}".format(round(multiplier,4),i+1,j+1))
                    showMatrix(C)  

    # Determining determinant

    det = 1

    for i in range(n):
        det = det*(C[i][i])

    for i in range(totalswaps):
        det = det*(-1)

    if(showall == 1):
        print("\nDeterminant of the Matrix A = {}\n".format(round(det,4)))

    eps = 0.0000001

    if(math.fabs(det) <= eps):
        print("Sorry, the given equations cannot be solved...\n")
        S = np.zeros((1,2))  # Will be diffentiated by shape
        return S

    # Performing Back substitution    
    
    if(showall == 1):
        print("Performing back substitution\n")

    S = np.zeros((n,1))  # Solution vector initialized

    for i in range(n-1,-1,-1):  # in reverse order
        p = C[i][n]
        for j in range(n-1,i,-1):
            p = p - C[i][j]*S[j][0]
        S[i][0] = p/C[i][i]

    
    # Returning solution vector

    return S

# returns m
# y = mx
def proportionalreg(xs,ys):
    
    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    sumxy = 0
    sumx2 = 0

    for i in range(n):
        sumxy = sumxy + xs[i]*ys[i]
        sumx2 = sumx2 + (xs[i]**2)

    m = sumxy/sumx2

# Linear Regression Start

# returns a and b
# Regression for y = a + bx
def linreg(xs,ys):

    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    sumxy = 0
    sumx = 0
    sumy = 0
    sumx2 = 0

    for i in range(n):
        sumxy = sumxy + xs[i]*ys[i]
        sumx = sumx + xs[i]
        sumy = sumy + ys[i]
        sumx2 = sumx2 + (xs[i]**2)

    b = (n*sumxy - (sumx*sumy))/(n*sumx2 - (sumx**2))
    a = (sumy/n) - b*(sumx/n)

    return a,b

def linreg_predict(xs,ys,allx):
    allx = np.array(allx)
    a,b = linreg(xs,ys)
    ally = a + b*allx
    return ally

# returns Sum of squares of residuals
def linreg_SSE(xs,ys):
    xs = np.array(xs)
    n = xs.shape[0]
    predicted_y = linreg_predict(xs,ys,xs)
    SSE = 0.0
    for i in range(n):
        SSE += (ys[i]-predicted_y[i])**2
    return SSE

def linreg_graph(xs,ys):  

    xs = np.array(xs)
    ys = np.array(ys)

    range = xs.max()-xs.min()

    allx = np.linspace(xs.min()-range/5,xs.max()+range/5,10000)
    ally = []

    plt.title("Plot for Regression")
    plt.plot(xs,ys,"og",label="Given Data Points")
    
    ally =  linreg_predict(xs,ys,allx)
    plt.plot(allx,ally,".-b",label="Predicted Curve",markersize=1)

    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.xlim([xs.min()-range/5,xs.max()+range/5])
    plt.legend(loc="best")
    plt.grid()
    plt.show()


# Linear Regression End

# Polynomial Regression Start

# y = polynomial of power m
# y = a0 + (a1)x + ..... + (am)x^m
# returns an array [a0 , a1 , a2 , a3 , ......... , am]
def polyreg(xs,ys,m):

    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    sumofxpowers = np.zeros(2*m+5)
    sumofxpowy = np.zeros(m+5)
    for i in range(n):
        x = xs[i]
        #print(x)
        tmp = 1.0
        for j in range(2*m+1):
            sumofxpowers[j] += tmp
            tmp = tmp*x
            #print(tmp)
        tmp = 1.0
        for j in range(m+1):
            sumofxpowy[j] += tmp*ys[i]
            tmp = tmp*x
    
    #print(sumofxpowers)
    #print(sumofxpowy)

    A = np.zeros((m+1,m+1))
    B = np.zeros((m+1,1))

    for i in range(m+1):
        for j in range(m+1):
            A[i][j] = sumofxpowers[i+j]
            B[i][0] = sumofxpowy[i]

    coeffs = GaussianElimination(A,B,1,0)

    #print(coeffs.transpose()[0])
    return coeffs.transpose()[0]

def polyreg_predict(xs,ys,allx,m):
    allx = np.array(allx)
    coeffs = polyreg(xs,ys,m)
    n = allx.shape[0]
    ally = np.zeros(n)
    for i in range(n):
        for j in range(m+1):
            ally[i] += coeffs[j]*(allx[i]**j)
    return ally

# returns Sum of squares of residuals
def polyreg_SSE(xs,ys,m):
    xs = np.array(xs)
    n = xs.shape[0]
    predicted_y = polyreg_predict(xs,ys,xs,m)
    SSE = 0.0
    for i in range(n):
        SSE += (ys[i]-predicted_y[i])**2
    return SSE

def polyreg_graph(xs,ys,m):  

    xs = np.array(xs)
    ys = np.array(ys)

    range = xs.max()-xs.min()

    allx = np.linspace(xs.min()-range/5,xs.max()+range/5,10000)
    ally = []

    plt.title("Plot for Regression")
    plt.plot(xs,ys,"og",label="Given Data Points")
    
    ally = polyreg_predict(xs,ys,allx,m)
    plt.plot(allx,ally,".-b",label="Predicted Curve",markersize=1)

    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.xlim([xs.min()-range/5,xs.max()+range/5])
    plt.legend(loc="best")
    plt.grid()
    plt.show()

# Polynomial Regression End

# Exponential Regression with Transformation Start

# y = a*e^(bx)
# returns a, b
# transforms to logarithmic scale, solves for linear regression and then to exponential
def expreg_transform(xs,ys):

    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    ys = np.log(ys)

    a,b = linreg(xs,ys)

    a = np.exp(a)

    return a,b

def expreg_transform_predict(xs,ys,allx):
    allx = np.array(allx)
    a,b = expreg_transform(xs,ys)
    n = allx.shape[0]
    ally = np.zeros(n)
    for i in range(n):
        ally[i] = a*np.exp(b*allx[i])
    return ally

# returns Sum of squares of residuals
def expreg_transform_SSE(xs,ys):
    xs = np.array(xs)
    n = xs.shape[0]
    predicted_y = expreg_transform_predict(xs,ys,xs)
    SSE = 0.0
    for i in range(n):
        SSE += (ys[i]-predicted_y[i])**2
    return SSE

def expreg_transform_graph(xs,ys):  

    xs = np.array(xs)
    ys = np.array(ys)

    range = xs.max()-xs.min()

    allx = np.linspace(xs.min()-range/5,xs.max()+range/5,10000)
    ally = []

    plt.title("Plot for Regression")
    plt.plot(xs,ys,"og",label="Given Data Points")
    
    ally = expreg_transform_predict(xs,ys,allx)
    plt.plot(allx,ally,".-b",label="Predicted Curve",markersize=1)

    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.xlim([xs.min()-range/5,xs.max()+range/5])
    plt.legend(loc="best")
    plt.grid()
    plt.show()
    
# Exponential Regression with Transformation End

# Exponential Regression Direct (without Transformation) Start

def expreg_direct(xs,ys):

    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    # guess for b
    loguess = -5 # subject to change
    higuess = 5 # subject to change
    err = 0.001
    mx = 10000

    def func(b):

        P = 0.0
        Q = 0.0
        R = 0.0
        S = 0.0

        for i in range(n):
            P += ys[i]*xs[i]*np.exp(b*xs[i])
            Q += ys[i]*np.exp(b*xs[i])
            R += np.exp(2.0*b*xs[i])
            S += xs[i]*np.exp(2.0*b*xs[i])
        
        return P - (Q*S)/R

    b = bisection(loguess,higuess,err,mx,func)
    Q = 0.0
    R = 0.0
    for i in range(n):
        Q += ys[i]*np.exp(b*xs[i])
        R += np.exp(2.0*b*xs[i])
    a = Q/R
    return a,b

def expreg_direct_predict(xs,ys,allx):
    allx = np.array(allx)
    a,b = expreg_direct(xs,ys)
    n = allx.shape[0]
    ally = np.zeros(n)
    for i in range(n):
        ally[i] = a*np.exp(b*allx[i])
    return ally

# returns Sum of squares of residuals
def expreg_direct_SSE(xs,ys):
    xs = np.array(xs)
    n = xs.shape[0]
    predicted_y = expreg_direct_predict(xs,ys,xs)
    SSE = 0.0
    for i in range(n):
        SSE += (ys[i]-predicted_y[i])**2
    return SSE

def expreg_direct_graph(xs,ys):  

    xs = np.array(xs)
    ys = np.array(ys)

    range = xs.max()-xs.min()

    allx = np.linspace(xs.min()-range/5,xs.max()+range/5,10000)
    ally = []

    plt.title("Plot for Regression")
    plt.plot(xs,ys,"og",label="Given Data Points")
    
    ally = expreg_direct_predict(xs,ys,allx)
    plt.plot(allx,ally,".-b",label="Predicted Curve",markersize=1)

    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.xlim([xs.min()-range/5,xs.max()+range/5])
    plt.legend(loc="best")
    plt.grid()
    plt.show()

# Exponential Regression Direct (without Transformation) End

# Tranformation to Any model from a given model

def anytransform(xs,ys):

    xs = np.array(xs)
    ys = np.array(ys)

    if(xs.shape[0] != ys.shape[0]):
        print("Sorry, Array dimensions not equal.")
        return
    
    n = xs.shape[0]

    newxs = np.zeros(n)
    newys = np.zeros(n)

    for i in range(n):
        newxs[i] = (np.exp(xs[i]))/xs[i]  # transform x
        newys[i] = ys[i]/xs[i]  # transform y

    print("Should be linear")
    print(newxs)
    print(newys)

    aprime , bprime = linreg(newxs,newys)  # Use tranformed regression

    a = aprime
    b = bprime

    return a,b

def anytransform_predict(xs,ys,allx):
    allx = np.array(allx)
    a,b = anytransform(xs,ys)
    n = allx.shape[0]
    ally = np.zeros(n)
    for i in range(n):
        ally[i] = a*allx[i] + b*np.exp(allx[i])  # write the given equation in terms of allx array
    return ally

# returns Sum of squares of residuals
def anytransform_SSE(xs,ys):
    xs = np.array(xs)
    n = xs.shape[0]
    predicted_y = anytransform_predict(xs,ys,xs)
    SSE = 0.0
    for i in range(n):
        SSE += (ys[i]-predicted_y[i])**2
    return SSE

def anytransform_graph(xs,ys,predictx):  

    xs = np.array(xs)
    ys = np.array(ys)

    range = xs.max()-xs.min()

    allx = np.linspace(xs.min()-range/5,xs.max()+range/5,10000)
    ally = []
    
    ally = anytransform_predict(xs,ys,allx)
    plt.title("Plot for Regression")

    plt.plot(allx,ally,".-b",label="Predicted Curve",markersize=1)

    plt.plot(xs,ys,"og",label="Given Data Points")

    predicty = anytransform_predict(xs,ys,predictx)

    plt.plot(predictx,predicty,"or",label="Predicted Points")

    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.xlim([xs.min()-range/5,xs.max()+range/5])
    plt.legend(loc="best")
    plt.grid()
    plt.show()



def compareModels(xs,ys):
    N = 80
    M = 5
    gap = " "
    spaces = 4*" " + "|" + 4*" "
    print("\n\nComparison among Models\n")
    print(N*"-")
    print(5*gap + "           Name of Model          " + spaces + gap + "Sum of Squared Residuals")
    print(N*"-")
    print(5*gap + "         Linear Regression        " + spaces + M*gap + "%.10f"%linreg_SSE(xs,ys))
    print(N*"-")
    print(5*gap + "   Polynomial Regression (m = 2)  " + spaces + M*gap + "%.10f"%polyreg_SSE(xs,ys,2))
    print(N*"-")
    print(5*gap + "   Polynomial Regression (m = 3)  " + spaces + M*gap + "%.10f"%polyreg_SSE(xs,ys,3))
    print(N*"-")
    print(5*gap + "Exponential Regression (transform)" + spaces + M*gap + "%.10f"%expreg_transform_SSE(xs,ys))
    print(N*"-")
    print(5*gap + "Exponential Regression ( direct ) " + spaces + M*gap + "%.10f"%expreg_direct_SSE(xs,ys))
    print(N*"-")
    print("\n")

def plot_data_points(xs,ys):  

    xs = np.array(xs)
    ys = np.array(ys)

    plt.title("Given Data Points")
    plt.plot(xs,ys,"o-g",label="Given Data Points")
    plt.xlabel('x-value')
    plt.ylabel('y-value')
    plt.legend(loc="best")
    plt.grid()
    plt.show()

#y = 0.11766514 + 0.09609143*x

"""
xvals = [0.698132, 0.959931, 1.134464, 1.570796, 1.919862]
yvals = [0.188224, 0.209138, 0.230052, 0.250965, 0.313707]
print(linreg(xvals,yvals))

"""

# y = 6.022e-6 + 6.278e-9x - 1.222e-11x^2

"""
xvals = [80,40,-40,-120,-200,-280,-340]
yvals = np.array([6.47,6.24,5.72,5.09,4.30,3.33,2.45])*(10**(-6))
print(polyreg(xvals,yvals,2))

"""

# y = 0.99974*e^(-0.11505x)

"""
xvals = [0,1,3,5,7,9]
yvals = [1,0.891,0.708,0.562,0.447,0.355]
print(expreg_transform(xvals,yvals))

"""

#xs = [0.5,0.8,1.5,2.5,4.0]
#ys = [1.1,2.4,5.3,7.6,8.9]

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


# Preparing Data for Numpy Float Array

xs = np.array(xs)
xs = np.asfarray(xs)

ys = np.array(ys)
ys = np.asfarray(ys)

#print(xs)
#print(ys)

# PLot given data points

plot_data_points(xs,ys)

#compareModels(xs,ys)
a,b = anytransform(xs,ys)
print("\nModel : " + str(a) + " x + " + str(b) + " e^x \n")


predictx = [1.45]

predicty = anytransform_predict(xs,ys,predictx)

for i in range(len(predictx)):
    print("Prediction for x = " + str(predictx[i]) + " :  y = " + str(predicty[i]))

print()

anytransform_graph(xs,ys,predictx)

# Always write the arrays in float





    


    
