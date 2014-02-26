import math

umin_sparse = 0.3
umax_sparse = 0.6
umin_common = 0.3
umax_common = 0.6
umin_pert = 0.3
umax_pert = 0.6
sparsity_prob = 0.98
percentage_common = 2
percentage_pert = 1

def m_common(p):
    math.ceil((percentage_common/100)*p)

def m_pert(p):
    math.ceil((percentage_pert/100)*p)
