import matplotlib.pyplot as plt
import numpy as np
import math

def f(x):
    a = 1
    b = -0.18
    c = 0
    d = 0.0004752
    return a*x*x*x + b*x*x + c*x + d

def fprime(x):   # The derivative of the function f(x)
    a = 1
    b = -0.18
    c = 0
    d = 0.0004752
    return 3*a*x*x + 2*b*x + c

def graph():
    x = np.linspace(-0.2,0.2,1000)
    y = f(x)
    plt.plot(x,y,"-g")
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.xlim([-0.15,0.2])
    plt.grid()
    plt.show()

def bisection(lo,hi,err,mx):
    print("\n\nBisection Method\n")
    print('\n     Iteration No.    |        x_l        |        x_m        |        x_u        |   Relative Approx Error   |        f(x_m)   ')
    ss = "-"
    for i in range(130):
        ss += '-'
    print(ss)
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
            if(cnt<10):
                ch1 = " "
            else:
                ch1 = ""
            if(e<10.0):
                ch = " "
            else:
                ch = ""
            if(fm>=0.0):
                ch2 = " "
            else:
                ch2 = ""
            print('        {}  {}          |     {:.7f}     |     {:.7f}     |     {:.7f}     |         {:.6f} %    {}   |    {}  {:.6f}     '.format(ch1,cnt,prelo,mid,prehi,e,ch,ch2,fm))
            #print("      " + str(cnt) +" Relative Approx. Error = {:.6f}%".format(e))
        else:
            print('           {}          |     {:.7f}     |     {:.7f}     |     {:.7f}     |            N/A            |       {:.6f}     '.format(cnt,prelo,mid,prehi,fm))
            #print("Iteration " + str(cnt) + ": Relative Approx. Error = Not Applicable")
        
        ans = mid

    ans = mid
    return ans

def newtonraphson(guess,err,mx):
    print("\nNewton-Raphson Method\n")
    print('\n     Iteration No.    |        x_i        |   Relative Approx Error   |      f(x_i)   ')
    offset = 0.001
    minrange = -10000
    maxrange = 10000
    ss = "-"
    for i in range(90):
        ss += '-'
    print(ss)
    eps = 10**(-18)
    cnt = 0
    e = 100 # dummy
    ans = guess
    while(e - err > (-eps)):
        if(cnt == mx):
            break
        while(math.fabs(fprime(ans) - 0.0) < 10**(-8)):  # Handling division by zero at minima and maxima
            ans  = ans + offset
        x = ans - (f(ans)/fprime(ans))

        # Handling divergence

        if(x > maxrange or x < minrange):
            print("\n\nSorry, Newton-Raphson is diverging. Please choose a new guess: ")
            newguess = int(input())
            newtonraphson(newguess,err,mx)
        
        e = (math.fabs((x-ans)/x))*100
        fx = f(x)
        cnt = cnt + 1
        if(cnt<10):
            ch1 = " "
        else:
            ch1 = ""
        if(e<10.0):
            ch = " "
        else:
            ch = ""
        if(fx>=0.0):
            ch2 = " "
        else:
            ch2 = ""
        print('          {}{}          |     {:.7f}     |         {:.6f} %   {}    |   {}   {:.6f}     '.format(ch1,cnt,x,e,ch,ch2,fx)) 
        #print("Iteration " + str(cnt) + ": Relative Approx. Error = {:.6f}%".format(e))
        ans = x
    return ans


graph()

low = 0.04
high = 0.08
guess = 0.05  # Upon observation of the graph
exp_err = 0.5
max_iter = 100

"""

print("Enter lower bound for bisection : ")
low = int(input())
print("Enter upper bound for bisection : ")
low = int(input())
print("Enter primary guess for Newton-Raphson : ")
low = int(input())

"""

ans = bisection( low , high , exp_err , max_iter )
print("\n\nSolution: x = " + str(ans) + "\n")

ans = newtonraphson( guess , exp_err , max_iter )
print("\n\nSolution: x = " + str(ans) + "\n\n")
