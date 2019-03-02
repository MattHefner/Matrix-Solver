

import math
global n, p, i, j, solutionvector, max, maxrow, temp, a1, a2, a3, x1, x2, x3

p = 0                 #We declare the variables that are going to be used and set them to a default value
i = 0                 # This allows the code to work for any size matrix that the user desires
j = 0
pivot = 0             # The variables p, i, and j are used when locating certain numbers in the matrix to do operations on
count = 0              # pivot is used to keep track of the pivot when row swapping and during operations
                       # count is a variable that is placed and increased everywhere a operation takes place and is necessary in counting operations
debug = False
# When debug = True, a sample matrix is inputed for testing purposes

def FindMatrix(A,n):               # This definition prompts the user to enter in specific data to build the matrix
    for i in range(n):
        for j in range(n+1):
            A[i][j] = int(input("Enter a number for row " + str(i+1) + ", column "+ str(j+1) + ":  "))
            
def RowChange(rows,m,n):         # RowChange is used when swapping the rows and is used in Part_Piv and Scale_Piv pivoting
    temp = rows[m]
    rows[m]=rows[n]
    rows[n] = temp 
              
def B_S(a,n, rows):           # B_S uses back substitution to solve for the solution vector 'x'
    global count
    x = [0 for i in range(n)]
    
    r = rows[n-1]
    x[n-1] = float(a[r][n])/a[r][n-1]              # B_S is calculating each value in the solution vector by solving the left side of the 
    count += 1                                     # equation and using operations to find each element. After finding the first value
                                                   # the code loops, plugs in the previous value where needed and solves until the whole vector is found
    for i in range (n-1,-1,-1):
      z = 0
      r = rows[i]
      for j in range(i+1,n):
         z = z  + float(a[r][j])*x[j]
         count += 2
         x[i] = float(a[r][n] - z)/a[r][i]
         count += 2
    return x

def Get_S(matrix,n,s):                         # Get_S finds the S values for Scale_Piv Part_Piv pivoting
    for i in range (0,n):
        max = 0
        for j in range (0, n):
            if max < math.fabs(matrix[i][j]):
                max = math.fabs(matrix[i][j])
                s[i] = max
    return s

def Gaus(a,p,n,rows,pivot):           # Gaus is a defintion that calculates the solved matrix through gaussian elimination
    global count
    for i in range(p+1,n):
        r = rows[i]                         
        m=-a[r][p]/a[pivot][p]          # We loop through the matrix and replace the values under the pivots with 0 and solve using operations
        count += 1
        a[r][p] = 0                         
        for j in range(p+1,n+1):            
          a[r][j] = float(a[r][j] + m * a[pivot][j])
          count += 2
              
def Part_Piv(a,p,n,rowVector, pivotRow):      # Part_Piv calculated the solved matrix by doing partial pivoting
    global count
    for p in range(0,n-1):                       
        max = 0                                
        for i in range(p,n):                  
            if abs(a[i][p]) > max:          # We are finding the max number in the columns to be the pivots and swap rows based on pivots
                max = abs(a[i][p])          # the defintion then calls in Gaus to solve the matrix
                maxrow = i                      
            else:
                maxrow = p

        if maxrow > p:                          
            RowChange(rowVector,p,maxrow)
            pivotRow = maxrow
        else:
            pivotRow = p
            
        Gaus(a,p,n,rowVector,pivotRow)

def Scale_Piv(a,p,n,rowVector,pivotRow,sVector):      #Scale_Piv is a defintion that does scaled partial pivoting
    global count                                      # This is similar to partial pivoting except that we know find the S vector, taking the
    s = Get_S(a,n,sVector)                             # max numbers from each row and using them later to divide by to find which rows need to swap
                                                       # Once the rows swap, gaussian elimination is used to find the solution matrix at each swap
    for p in range(0,n-1):                      
        max = 0                                 
        for i in range(p,n):                    
            row = rowVector[i]
            m[i] = math.fabs(a[row][p])/s[row]
            count += 1
            if m[i] > max:
                max = m[i]
                maxrow = i
                
        if maxrow > p:                               # Evaluates and finds where the pivot should be and executes RowChange to swap
            RowChange(rowVector,p,maxrow)
            pivotRow = maxrow
        else:
            pivotRow = rowVector[p]       
            
        Gaus(a,p,n,rowVector,pivotRow)


selectedMethod = int(input("Select a method to solve the system of equations (1: Gaussian, 2: Partial Pivoting, 3: Scaled Partial Pivoting): "))

n = int(input("Enter size of a square matrix (optimal results for sizes >= 3): "))          #This selected method calls Gaus to solve the matrix, the user is also
a = [[0 for x in range(n+1)] for y in range(n)]                   # prompted to enter in the size matrix to cut back on errors
s = [0 for x in range(n)]
m = [0 for x in range(n)]
Part_PivRV = [i for i in range(n)]

if debug == True:
    a1 = [2,3,1,-4]
    a2 = [4,1,4,9]            # this is the debug matrix used for testing
    a3 = [3,4,6,0]
    a = [a1, a2, a3]
else:
    FindMatrix(a,n)



if selectedMethod == 1:  
    for p in range (0,n):
        Gaus(a,p,n,Part_PivRV, Part_PivRV[p])        
    
    solutionvector = B_S(a,n,Part_PivRV)
        
    print("\n\nSolved Matrix:\n" + str(a[0]) + "\n" + str(a[1]) + "\n" + str(a[2]))            #This outputs the information to the users
    print("\nSolution Vector: " + str(solutionvector) + ", with " + str(count) + " function evaluations.")
    
if selectedMethod == 2:
                               #Calls Part_Piv to solve the matrix and outputs the answers 
    
    Part_Piv(a,p,n,Part_PivRV, pivot)
        
    solutionvector = B_S(a,n,Part_PivRV)
    
    print("\n\nSolved Matrix:\n" + str(a[0][0:n+1]) + "\n" + str(a[1][0:n+1]) + "\n" + str(a[2][0:n+1]))
    print("\nSolution Vector: " + str(solutionvector) + ", with " + str(count) + " function evaluations.")
    
if selectedMethod == 3:
                               #Calls Scale_Piv to solve the matrix and outputs the answers
    Scale_Piv(a,p,n,Part_PivRV,pivot,s)
    
    solutionvector = B_S(a,n,Part_PivRV)
    
    print("\n\nSolved Matrix:\n" + str(a[0][0:n+1]) + "\n" + str(a[1][0:n+1]) + "\n" + str(a[2][0:n+1]))
    print("\nSolution Vector: " + str(solutionvector) + ", with " + str(count) + " function evaluations.")