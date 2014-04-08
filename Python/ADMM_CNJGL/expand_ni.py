from numpy.linalg import eigh
from numpy import diag, transpose, sqrt

# B: square matrix
# S: square matrix

def expand_ni(B, S, rho, n):
    A = B - S*n / (2.0 * rho)

    (S, U) = eigh(A, 'U')

    s = S
    s = 0.5 * (s + sqrt( s*s + 2.0/(rho/(n+0.0)) ))

    mul1 = U*diag(s)
    ret = mul1 * transpose(U)
    return ret

''' # test:
from numpy import matrix
print ( expand_ni(matrix([[1,2,3],[4,5,6],[7,8,9]]), 
                  matrix([[2,3,4],[5,2,4],[1,3,5]]), 2.0, 3.0) )
'''
