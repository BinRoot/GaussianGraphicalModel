	
% --------------------- LAST UPDATE: 8/1/2013 ------------------------------- %

function [Theta_true_1, Theta_true_2, ind_m_common, ind_m_pert] = data(p)

data_parameters; % Load data parameters 
flip = true;
%-------------------------------------------------------------------------
% Data Generation Steps: 

% [1] Build the Theta_true_1 and Theta_true_2 matrices, in 3 Layers
%       (1) Layer 1:  A sparse Matrix 
%       (2) Layer 2:  perturbed nodes 
%       (3) Layer 3:  common co-hub nodes

% [2] Conditioning: The data must have the following Characteristics: 
%       1. Symmetric
%       2. Positive Semi Definite 
%-------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Layer 1: A sparse Matrix 

% 1. Make a dense matrix with values drawn from U[umin_sparse, umax_sparse] 
% idea: make a base of umin_sparsevalues, then add any number, 
% between the values of umax_sparse and umax_sparse
dense = (ones(p,p) * umin_sparse) + (rand(p,p) * (umax_sparse - umin_sparse));

% 2. Make the "dense matrix" sparse. 
indmat = (rand(p,p) < sparsity_prob);
dense(indmat) = 0; 

Net_1 = dense;
Net_2 = dense;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (2) Layer 2:  Put in perturbed nodes into Net_1  

ind_m_pert = 0;
if m_pert > 0
 
    % randomly select "m_pert" number of cols to be perturbed in Net_1
ind_m_pert = unidrnd(p, 1, m_pert);

for i = ind_m_pert;
    for j = 1:p;
        Net_1(j, i) = umin_pert + (rand(1) * (umax_pert - umin_pert));
               
        Net_1(i, j) = umin_pert + (rand(1) * (umax_pert - umin_pert));
    end 
end 

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Layer 3:  Put in common co-hub nodes   (identical in Net_1 and Net_2)

% select "m_common" number of cols to be made common in both Net_1 and Net_1



ind_m_common = unidrnd(p, 1, m_common);

for i = ind_m_common;
    for j = 1:p;
        a = umin_common + (rand(1) * (umax_common - umin_common)); 
        Net_1(i,j) = a;
        Net_2(i,j) = a;
        
        Net_1(j,i) = a;
        Net_2(j,i) = a;
    end 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final Step, make all distributions "U[-umin, -umax] union U[umin, umax]",
% include negative parts. (No reason to do this sooner.) 

% Randomly select, 50% of values, to be negative. 
if (flip)
    flipper = unidrnd(2,p,p);
    flipper(flipper == 2) = -1;
    Net_1 = Net_1 .* flipper;
    Net_2 = Net_2 .* flipper;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [2] "Conditioning" the data; it must have the following Characteristics: 

% 1. Symmetric
Net_1_sym = triu(Net_1) + triu(Net_1)'; 
Net_2_sym = triu(Net_2) + triu(Net_2)';

    % check
issym=@(x) all(all(tril(x)==triu(x).'));
if ~issym(Net_1_sym)
    error('input Theta_1 is not symmetric.');
end
if ~issym(Net_2_sym)
    error('input Theta_2 is not symmetric.');
end


% 2. Positive Semi Definite 
min_eigval_1 = min(eig(Net_1_sym));
min_eigval_2 = min(eig(Net_2_sym));
    
min_eigval = min(min_eigval_1, min_eigval_2); 

    % adjust the diagonals 
Theta_true_1 = Net_1_sym + (eye(p)*(abs(min_eigval) + 0.1)); 
Theta_true_2 = Net_2_sym + (eye(p)*(abs(min_eigval) + 0.1)); 

    % check
if (min(eig(Theta_true_1)) <= 0)
    error('input Theta_1 is not PSD.');
end
if (min(eig(Theta_true_2)) <= 0)
    error('input Theta_2 is not PSD.');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end 
