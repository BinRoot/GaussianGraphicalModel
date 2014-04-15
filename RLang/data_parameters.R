umin_sparse <- 0.3 # Lower and upper bounds for magnitude of entries for the sparse matrix
umax_sparse <- 0.6
umin_common <- 0.3 # Lower and upper bounds for magnitude of entries in common cohub nodes
umax_common <- 0.6
umin_pert <- 0.3 # Lower and upper bounds for magnitude of perturbations
umax_pert <- 0.6
sparsity_prob <- 0.98 # Sparsity level desired in the true inverse covariance 
percentage_common <- 2
percentage_pert <- 1
m_common <- ceiling((percentage_common/100)*p)
m_pert <- ceiling((percentage_pert/100)*p)
