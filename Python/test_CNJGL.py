# Test file for ADMM algorithm

from time import time

import sys
sys.path.append('./ADMM_CNJGL')

from numpy import matrix, cov
from ADM_CNJGL import ADM_CNJGL
import problem_parameters

# Input data
X = matrix([[0.05, -0.04, 0.06, 0.02, -0.01],
            [0.02, 0.01, -0.06, 0.01, -0.04],
            [0.05, 0.04, -0.02, 0.01, -0.01]])

(xrows, xcols) = X.shape

S1 = cov(X)

S2 = S1


tic = time()

ADM_CNJGL(S1, S2, problem_parameters.lambda_1, problem_parameters.lambda_2, problem_parameters.n1, problem_parameters.n2)

elapsed = time() - tic

print "runtime: ", elapsed
