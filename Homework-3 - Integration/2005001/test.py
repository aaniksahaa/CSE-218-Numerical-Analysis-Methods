import numpy as np
import matplotlib.pyplot as plt
import math

cme = 5e-4

def func(x):
    up = 6.73*x + 6.725e-8 + (7.26e-4)*cme
    down = (3.62e-12)*x + (3.908e-8)*x*cme
    return (-up)/down

# n-segment trapezoidal rule integration
# argument f is the function to integrate
def trapezoid_integrate(a,b,n,f): # f is the integration function
    h = (b-a)/n
    ans = f(a)
    for i in range(1,n):
        ans = ans + 2*f(a+i*h)
    ans = ans + f(b)
    ans = ((b-a)*ans)/(2*n)
    return ans

# Shows the integral values and relative approximate error for
# Multiple Segment Trapezoidal Integration
# argument f is the function to integrate
def trapezoid(a,b,maxn,f):
    gap = " "
    spaces = 4*" " + "|" + 4*" "
    print("\n\nMultiple Segment Trapezoidal Integration\n")
    print(100*"-")
    print(5*gap + "Number of Segments" + spaces + gap + "Integral Value" + spaces + "Relative Approx. Error")
    print(100*"-")
    prev = -1 # dummy
    for i in range(1,maxn+1):
        answer = trapezoid_integrate(a,b,i,f)
        if(i>1):
            err = math.fabs((answer - prev)/answer)*100
        firstpart = 14*gap + str(i) + 8*gap + spaces + "%.6f"%answer + spaces + 4*gap
        if(i>1):
            print(firstpart + "%.10f"%err + " %")
        else:
            print(firstpart + "   -------   ")
        prev = answer

# n = number of segments
# divides the interval into n equal segments
# n must be even
def simpson_integrate(a,b,n,f):
    h = (b-a)/n
    ans = f(a)
    for i in range(1,n):
        if((i%2)==1):
            ans += 4*(f(a+i*h))
        else:
            ans += 2*(f(a+i*h))
    ans  = ans + f(b)
    ans = ((b-a)*ans)/(3*n)
    return ans

# maxn = maximum number of applications
def simpson(a,b,maxn,f):
    gap = " "
    spaces = 4*" " + "|" + 4*" "
    print("\n\nMultiple Application Simpson's 1/3 rule Integration\n")
    print(100*"-")
    print(5*gap + "Number of Segments" + spaces + gap + "Integral Value" + spaces + "Relative Approx. Error")
    print(100*"-")
    prev = -1 # dummy
    for i in range(2,maxn+1,2):
        answer = simpson_integrate(a,b,i,f)
        if(i>2):
            err = math.fabs((answer - prev)/answer)*100
        firstpart = 14*gap + str(i) + 8*gap + spaces + "%.6f"%answer + spaces + 4*gap
        if(i>2):
            print(firstpart + "%.10f"%err + " %")
        else:
            print(firstpart + "   -------   ")
        prev = answer

# n = number of segment
# method = trapezoid_integrate or simpson_integrate
def graph(initial,xvals,f,method,n):
    xvals = np.array(xvals)
    yy = []
    for x in xvals:
        yy.append(method(initial,x,n,f))
    yvals = np.array(yy)
    plt.plot(yvals,xvals,"-og")
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Oxygen Concentration")
    plt.show()


initial_oxygen = 1.22e-4

primary_fraction = 1
final_fraction = 0.5

x1 = primary_fraction*initial_oxygen
x2 = final_fraction*initial_oxygen

def g(x):
    return 2000*np.log(140000/(140000-2100*x))-9.8*x

x1 = 8
x2 = 30
n = 10

trapezoid(x1,x2,n,g)
simpson(x1,x2,n,g)

