# Test file for ADMM algorithm

from time import time
from numpy import matrix, transpose
import sys
sys.path.append('./ADMM_CNJGL')

from numpy import matrix, cov
from ADM_CNJGL import ADM_CNJGL
import problem_parameters

# Input data
X = matrix([[0.05, -0.04, 0.06, 0.02, -0.01],
            [0.02, 0.01, -0.06, 0.01, -0.04],
            [0.05, 0.04, -0.02, 0.01, -0.01]])

'''
import scipy.io
X = scipy.io.loadmat('./data/X.mat')
X = matrix(X['X'])[:,:500]
print X.shape
'''

(xrows, xcols) = X.shape

S1 = cov(transpose(X))
print "S1", S1.shape

S2 = S1


tic = time()

(Theta1, Theta2, V1, V2, iter_adm, relError) = ADM_CNJGL(S1, S2, problem_parameters.lambda_1, problem_parameters.lambda_2, problem_parameters.n1, problem_parameters.n2)

elapsed = time() - tic

print "runtime: ", elapsed

print Theta1
