# ADMM algorithm for CNJGL

import ADM_parameters
from numpy.linalg import norm
from expand_ni import expand_ni
from soft_thresh import soft_thresh
from soft_scal import soft_scal
from diag_construction import diag_construction
from nodiag_construction import nodiag_construction

from numpy import eye, zeros, diag, matrix, concatenate

def ADM_CNJGL(S1, S2, lam1, lam2, n1, n2):
    # Begin algorithm
    p = S1.shape[0]
    
    # Initialize variables
    Theta1old = eye(p)
    Theta2old = eye(p)
    Theta1 = eye(p)
    Theta2 = eye(p)
    V1 = Theta1old/2
    V2 = Theta2old/2
    W1 = V1.transpose()
    W2 = V2.transpose()
    Z1 = Theta1
    Z2 = Theta2
    Gamma1 = matrix(zeros((p,p)))
    Gamma2 = matrix(zeros((p,p)))
    Lambda1 = matrix(zeros((p,p)))
    Lambda2 = matrix(zeros((p,p)))
    Q1 = matrix(zeros((p,p)))
    Q2 = Q1

    # Load parameters
    rho = ADM_parameters.opts['rho_init']
    iter_adm = matrix(zeros((ADM_parameters.opts['homotopy_size'], 1)))
    
    relError = matrix(zeros((ADM_parameters.opts['maxiter'], 1)))

    # Begin outer loop for homotopy parameter
    for r in range(0, int(ADM_parameters.opts['homotopy_size']-1)):
        # Begin ADM
        for i in range(0, ADM_parameters.opts['maxiter'] - 1):
            # Update of Theta1 and Theta2
            Theta1 = expand_ni(1/(2*rho)*(rho*(V1 + W1 + Z1) - (Gamma1 + Lambda1)), S1, rho, n1)
            Theta2 = expand_ni(1/(2*rho)*(rho*(V2 + W2 + Z2) - (Gamma2 + Lambda2)), S2, rho, n2)
            
            # Update of Z1 and Z2
            Z1 = soft_thresh(Theta1 + Lambda1/rho, lam1/rho)
            Z2 = soft_thresh(Theta2 + Lambda2/rho, lam1/rho) # lam2?
            
            # Update of V1 and V2
            C1 = 1/(2*rho)*(rho*(W1.transpose() + Theta1 - W1) - (Q1 - Gamma1))
            C2 = 1/(2*rho)*(rho*(W2.transpose() + Theta2 - W2) - (Q2 - Gamma2));
            d1 = diag(C1)
            d2 = diag(C2)
            N1 = nodiag_construction(C1)
            N2 = nodiag_construction(C2)
            N = concatenate((N1,N2))
            H = soft_scal(N,lam2/(2*rho))
            d = concatenate((d1,d2))
            print "not dead 111"
            V1 = diag_construction(H[0:p-1,:],d1);
            print "not dead 222"
            V2 = diag_construction(H[p:2*p-3,:],d2); # err?
            print "not dead 333"
            
            # Update of W1 and W2
            W1 = 1/(2*rho)*(rho*(V1.transpose() + Theta1 - V1) + Gamma1 + Q1.transpose());
            W2 = 1/(2*rho)*(rho*(V2.transpose() + Theta2 - V2) + Gamma2 + Q2.transpose());

            # Update of dual variables
            Gamma1 = Gamma1 + rho*(Theta1 - V1 - W1);
            Gamma2 = Gamma2 + rho*(Theta2 - V2 - W2);
            Lambda1 = Lambda1 + rho*(Theta1 - Z1);
            Lambda2 = Lambda2 + rho*(Theta2 - Z2);
            Q1 = Q1 + rho*(V1 - W1.transpose());
            Q2 = Q2 + rho*(V2 - W2.transpose());

            relError[i + r*ADM_parameters.opts['maxiter'], 1] = \
                max(norm(Theta1 - Theta1old, 'fro') / norm(Theta1, 'fro'), \
                    norm(Theta2 - Theta2old, 'fro') / norm(Theta2, 'fro'))
            
            if relError[i + (r*ADM_parameters.opts['maxiter']), 1] <= ADM_parameters.opts['eps']:
                break
            
            Theta1old = Theta1
            Theta2old = Theta2
            # End ADM

    rho = rho * ADM_parameters.opts['homotopy']
    iter_adm[r,1] = i
    return (Theta1, Theta2, V1, V2, iter_adm, relError)
# End algorithm
