# Soft-scaling operator with diagonals
# Not penalized

from numpy import zeros, diag, matrix
from soft_scal import soft_scal
from diag_construction import diag_construction
from nodiag_construction import nodiag_construction

def soft_scal_nodiag(Y, lam):
    p = Y.shape[1]
    X = matrix(zeros((p,p)))
    xtilde = matrix(zeros((p-1,p-1)))
    ytilde = nodiag_construction(Y)
    
    xtilde = soft_scal(ytilde, lam)
    return diag_construction(xtilde, diag(Y))

