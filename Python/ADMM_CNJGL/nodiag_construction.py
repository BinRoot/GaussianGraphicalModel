# Remove diagonal elements of a matrix and rearrange
# A p x p matrix X is taken as input, its diagonals are removed.
# The resulting matrix is of size p-1 x p

from numpy import zeros, concatenate, matrix

def nodiag_construction(X):
    p = X.shape[0]
    xtilde = matrix(zeros((p-1,p-1)))
    
    xtilde[:,0] = X[1:p,0]

    for j in range(1, p-1):
        xtilde[:,j] = concatenate([X[0:j,j], X[j+1:p,j]])

    xtilde[:,p-2] = X[0:p-1, p-1]
    return xtilde

''' # test:
from numpy import matrix
A = matrix([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
print(nodiag_construction(A))
'''
