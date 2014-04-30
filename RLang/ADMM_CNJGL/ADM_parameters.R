# SET PARAMETERS FOR ADM ALGORITHM

opts.rho_init <- 1 # Initial rho
opts.maxiter <- 500 # Max iter
opts.eps <- 1e-6 # Termination tolerance
opts.homotopy <- 5 # Scaling parameter for rho
opts.rho_max <- 625

if (opts.rho_max/opts.rho_init == 1) {
    opts.homotopy_size <- 1
} else {
    opts.homotopy_size <- ceiling(log(opts.rho_max)/log(opts.homotopy)) + 1
}


