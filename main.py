import numpy as np
import math as mt
from functools import reduce

def find_v(A, step) :
    v = np.empty([len(A), 1])
    C = 0
    for i in range(0, step):
        A[i][step] = 0
    for i in range(len(A)):
       if i == step:
           tmp = reduce(lambda x, y: x+y, map(lambda x: x**2, np.transpose(A)[i]))
           C = np.sign(A[i][i])
           if C == 0:
               C = 1
           v[i] = A[i][i] + C * np.sqrt(tmp)
       elif i < step:
           v[i] = 0
       else:
           v[i] = np.transpose(A)[step][i]
    return v

def Diagonal(R):
    for i in range(len(R)):
        for j in range(len(R)):
            if j < i and np.fabs(R[i][j]) >= 0.01:
                return False
    return True

def QR(A):
    step = 0
    E = np.eye(len(A))
    A_s = A
    Q = np.eye(len(A))
    R = A
    while step < len(A) - 1:#not Diagonal(R):
        v = find_v(R, step)
        H = E - (2 / (np.dot(np.transpose(v), v))) * np.dot(v, np.transpose(v))
        R = np.dot(H, A_s)
        A_s = np.copy(R)
        Q = np.dot(Q, H)
        step += 1
    return (Q, R)

def AllFound(selfNumbers):
    for i in range(len(selfNumbers)):
        if selfNumbers[i] == 0 :
            return False
    return True

def ZeroColumn(A, i, eps) :
    column_eps = 0
    for j in range(i+1, len(A)):
        column_eps += A[j][i] * A[j][i]
    column_eps = mt.sqrt(column_eps)
    if column_eps <= eps :
        return True
    else :
        return False

def ComplexEmpty(compexPairs, i):
    complexPairs[i] = -1
    complexPairs[i + 1] = -1


def isComplexConcat(complexPairs, i):
    if complexPairs[i] == -1 and complexPairs[i+1] == -1:
        return True
    elif complexPairs[i] == i+1 and complexPairs[i+1] == i:
        return True
    else:
        return False

def MaybeComplex(A, i, complexPairs, eps,it) :
    #if complexPairs[i] == -1:
     #   complexPairs[i] = i+1
      #  complexPairs[i+1] = i
    if isComplexConcat(complexPairs, i):
        j = i + 2
        if j < len(A) :
            cur_eps = 0
            for x in range(j, len(A)) :
                cur_eps += A[x][i] * A[x][i]
            cur_eps = mt.sqrt(cur_eps)
            if(cur_eps <= eps and mt.fabs(A[i+1][i]) > eps) :
                complexPairs[i] = i + 1
                complexPairs[i+1] = i
                return True
            else:
                ComplexEmpty(complexPairs, i)
                return False
        elif j == len(A):
            if(mt.fabs(A[i+1][i]) > eps) :
                complexPairs[i] = i + 1
                complexPairs[i + 1] = i
                return True
            else:
                ComplexEmpty(complexPairs, i)
                return False
        else:
            ComplexEmpty(complexPairs, i)
            return False
    else:
        return False
def ComplexRoots(A, i) :
    a = 1
    b = -(A[i][i] + A[i+1][i+1])
    c = A[i][i]*A[i+1][i+1] - A[i][i+1]*A[i+1][i]
    roots = np.roots([a, b, c])
    #if isinstance(roots[0], complex):
        #return (roots[0], roots[1])
    #else :
        #return (0, 0)
    return (roots[0], roots[1])

def SelfFound(A, selfNumbers, complexDifference, complexPairs, eps, it):
    nextIter = False
    for i in range(len(A)) :
        if nextIter == True:
            nextIter = False
        else:
            if ZeroColumn(A, i, eps) == True and complexPairs[i] == -1:#not isinstance(selfNumbers[i], complex):
                selfNumbers[i] = A[i][i]
                complexDifference[i] = 0
            else:
                if MaybeComplex(A, i, complexPairs, eps, it) == True:
                    (l1, l2) = ComplexRoots(A, i)
                    complexDifference[i] = abs(selfNumbers[i] - l1)
                    complexDifference[i + 1] = abs(selfNumbers[i + 1] - l2)
                    selfNumbers[i] = l1
                    selfNumbers[i + 1] = l2
                    nextIter = True
                    #this was removing now upper is testing
                    #if isinstance(l1, complex):
                     #   nextIter = True
                      #  complexDifference[i] = abs(selfNumbers[i] - l1)          #new
                       # complexDifference[i + 1] = abs(selfNumbers[i + 1] - l2)  #new
                    #else :
                     #   complexPairs[i] = -1
                      #  complexPairs[i+1] = -1
                       # complexDifference[i] = eps+1                                    #new
                        #complexDifference[i + 1] = eps+1                                #new

                    #selfNumbers[i] = l1                                                 #new
                    #selfNumbers[i + 1] = l2                                             #new

def SmallDifference(complexDifference, eps):
    for i in range(len(complexDifference)):
        if complexDifference[i] > eps:
            return False
    return True



A = np.loadtxt("input.txt")
#A = np.array([[1, 2, 2, 8, 0],[1, 1, 5, 5, 2],[1, 51, 5, 11, 5], [1, 0, 0, 1, 3], [5, 0, 1000, 1, 212]])
Q = np.eye(len(A))
R = A
eps = 0.01
selfNumbers = [0] * len(A)
complexDifference = [eps+1] * len(A)
complexPairs = [-1] * len(A)
it = 0
while not AllFound(selfNumbers) or not SmallDifference(complexDifference, eps):
    (Q, R) = QR(A)
    A = np.dot(R, Q)
    SelfFound(A, selfNumbers, complexDifference, complexPairs, eps, it)
    it += 1
selfNumbers = np.round(selfNumbers, decimals=2)
print(selfNumbers)
print(np.round(np.linalg.eigvals(A),decimals=2))
np.savetxt("out.txt", selfNumbers, fmt="%.2f   ")
