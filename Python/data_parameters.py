from math import ceil

# Lower and upper bounds for magnitude of entires for the sparse matrix
umin_sparse = 0.3
umax_sparse = 0.6

# Lower and upper bounds for magnitude of entires in common cohub nodes
umin_common = 0.3
umax_common = 0.6

# Lower and upper bounds for magnitude of perturbations
umin_pert = 0.3
umax_pert = 0.6

# Sparsity level desired in the true inverse covariance matrix
sparsity_prob = 0.98


percentage_common = 2
percentage_pert = 1

# Number of common cohub nodes
def m_common(p):
    ceil((percentage_common/100)*p)

# Number of perturbed nodes
def m_pert(p):
    ceil((percentage_pert/100)*p)
