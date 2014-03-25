from math import sqrt
from numpy.linalg import eig
from numpy import diag, transpose

# B: square matrix
# S: square matrix

def expand_ni(B, S, rho, n):
    A = B - S*n / (2.0 * rho)
    (S, U) = eig(A)
    s = 0.5 * (S + map( lambda x: sqrt(x), S*S + 2/(rho/n) ))
    return U * diag(s) * transpose(U)

''' # test:
from numpy import matrix
print ( expand_ni(matrix([[1,2,3],[4,5,6],[7,8,9]]), 
                  matrix([[2,3,4],[5,2,4],[1,3,5]]), 2.0, 3.0) )
'''
