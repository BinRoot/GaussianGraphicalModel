%% ------------------- SOFT-SCALING OPERATOR WITH DIAGONALS   ---------- %%
%% ------------------- NOT PENALIZED ----------------------------------- %%

% -------------------- KARTHIK MOHAN, EE, UW -------------------------%
% -------------------- LAST UPDATE: 8/15/2012 ------------------------%

% ---------------- minimize_X \frac{1}{2}||X - Y||_2^2 
% ----------------             + \lambda*sum_i||\tilde{X}_i||_2
% ---------------- where \tilde{X} is X with diagonals removed and
% ---------------- rearranged into a square matrix.

% ---------------- solution is  \tilde{X}_i = 
% ---------------- max(||\tilde{Y}_i||_2 - \lambda, 0)*\tilde{Y}_i/||\tilde{Y}_i||_2
% ---------------- Also, X_{ii} = Y_{ii} \forall i.



function [X] = soft_scal_nodiag(Y,lambda)

p = size(Y,2);
X = zeros(p,p);
Xtilde = zeros(p-1,p-1);
Ytilde = nodiag_construction(Y);

Xtilde = soft_scal(Ytilde,lambda);
X = diag_construction(Xtilde,diag(Y));



