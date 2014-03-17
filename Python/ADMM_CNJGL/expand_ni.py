import numpy

def expand_ni(B, S, rho, n):
    A = B - S*n / (2.0 * rho)
    (S, U) = numpy.linalg.eig(A)
    s = 0.5 * (S + sqrt(S*S + 2/(rho/n)))
    return U * numpy.diag(s) * numpy.transpose(U)
    
