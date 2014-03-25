# Soft-scaling operator

from numpy.linalg import norm
from numpy import zeros, matrix

def soft_scal(Y, lam):
    m = Y.shape[0]
    n = Y.shape[1]
    X = matrix(zeros((m,n)))
    for t in range(0, n):
        if norm(Y[:,t]) > lam:
            X[:,t] = norm(Y[:,t] - lam) / norm(Y[:,t]) * Y[:,t].getA1()
        else:
            X[:,t] = matrix(zeros((m,1)))
    return X


''' # test:
import numpy
myY = numpy.matrix([[1,2,3],[4,5,6],[7,8,9]])
myLam = 3
print ( soft_scal(myY, myLam) )
'''
