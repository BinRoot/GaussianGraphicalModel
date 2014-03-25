# Soft-thresholding operator

from math import copysign

def soft_thresh(Y, lam):

    signY = Y
    # perform sign on matrix
    for c in range(0, Y.shape[0]):
        for r in range(0, Y.shape[1]):
            if Y[r,c] != 0:
                signY[r,c] = copysign(1, Y[r,c])

    AbsYminusLam = abs(Y) - lam
    for c in range(0, AbsYminusLam.shape[0]):
        for r in range(0, AbsYminusLam.shape[1]):
            if AbsYminusLam[r,c] < 0:
                AbsYminusLam[r,c] = 0

    return signY * AbsYminusLam
    
''' # test:
from numpy import matrix
myY = matrix([[1,2],[0,-4]])
print(soft_thresh(myY, 0))
'''
