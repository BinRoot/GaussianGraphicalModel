%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------- SYNTHETIC DATA GENERATION PARAMETERS ------------------%
%									     %	
% LAST UPDATE: 3/1/2013 						     %	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ----------- DATA generation parameters
umin_sparse = 0.3; % Lower and upper bounds for magnitude of entries for the sparse matrix
umax_sparse = 0.6;
umin_pert = 0.3; % Lower and upper bounds for magnitude of perturbations
umax_pert = 0.6;
sparsity_prob = 0.9; % Sparsity level desired in the true inverse covariance matrix
 

percentage_perturbed = 2;
m_pert = ceil((percentage_perturbed/100)*p); % Number of perturbed nodes

