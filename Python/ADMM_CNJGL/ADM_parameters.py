from math import ceil, log

opts = {
    'rho_init': 1, # Initial rho
    'maxiter': 500, # Max iter
    'eps': 1e-6, # Termination tolerance
    'homotopy': 5, # Scaling parameter for rho
    'rho_max': 625,
}

if opts['rho_max'] / opts['rho_init'] == 1 :
    opts['homotopy_size'] = 1
else :
    opts['homotopy_size'] = ceil ( log(opts['rho_max']) / 
                                   log(opts['homotopy']) ) + 1
