import math

opts = {
    'rho_init': 1,
    'maxiter': 500,
    'eps': 1e-6,
    'homotopy': 5,
    'rho_max': 625,
}

if opts['rho_max'] / opts['rho_init'] == 1 :
    opts['homotopy_size'] = 1
else :
    opts['homotopy_size'] = math.ceil ( math.log(opts['rho_max']) / 
                                        math.log(opts['homotopy'])) + 1
