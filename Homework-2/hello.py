import numpy as np

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

def  GaussianElimination(A,B,pivot,showall):

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
                    print("   First element is zero. Finding non-zero element...")
                
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

    if(det == 0):
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

print("\nEnter the number of unknown variables : \n")
n = int(input())

A = []

print("\nInput the matrix A : \n")

for i in range(n):
    row = list(map(float, input().split()))
    A.append(row)

print("\nInput the matrix B :\n")

B = []

for i in range(n):
    temp = []
    p = float(input())
    temp.append(p)
    B.append(temp)

showoptions = 1

# default

pivot = 1
showall = 1

# Input the choice for pivot and showall

if(showoptions == 1):

    print("\nDo you want partial pivoting? (Type 0 or 1)\n")
    pivot = int(input())

    print("\nDo you want demonstration of steps? (Type 0 or 1)\n")
    showall = int(input())

print()
print()

# Perform Gaussian Elimination

solution_vector = GaussianElimination(A,B,pivot,showall)

# Display the solutions

if(solution_vector.shape[1] == 1):
    print("Solution to the equations (Upto 4 decimal places): \n")
    for i in range(n):
        for j in range(1):
            print("{:.4f}".format(solution_vector[i][j]),end="    ")
        print()
    print()