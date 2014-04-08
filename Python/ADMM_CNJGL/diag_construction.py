# Takes a p-1 x p matrix and appends diagonals to get square matrix
#   Let A be p-1 x p matrix. Let d be a vector of p entires.
#   These d elements are added to the diagonal positions
#   of A to create a p x p matrix.

from numpy import zeros, matrix

# xtilde is a numpy matrix, d is an array
def diag_construction(xtilde, d):
    p = len(d)
    X = matrix(zeros((p, p)))

    X[0,0] = d[0]
    X[1:p,0] =  xtilde[:,0]

    for j in range (1, p-1):
        X[0:j,j] = xtilde[0:j,j]
        X[j,j] = d[j]
        X[j+1:p,j] = xtilde[j:p,j]

    X[0:p-1,p-1] = xtilde[:,p-1]
    X[p-1,p-1] = d[p-1]
    return X

''' # test: 
from numpy import matrix
xt = matrix([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
d = [99,99,99,99]
v = diag_construction(xt,d)
print(v)
'''
