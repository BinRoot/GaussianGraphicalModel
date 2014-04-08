%% ----------- The Expand Operator ----------------------- %%

% --- Expand(B, rho) = argmin_{Theta} n(-\log\det(Theta) + <Theta,S>) + rho*||Theta -
% B||_F^2

% --- Let A = B - S*n/(2*rho) = U*Lambda*U' be the Eigen-Decomposition of A.
% Then Set Lambda_{ii} = 0.5*(Lambda_{ii} + sqrt(Lambda_{ii}^2 + 2/rho));
% Expand(A, rho) = U*Sigma2*V'.




function [Theta] = expand_ni(B,S,rho,n)
A = B - S*n/(2*rho);

[U,S] = eig(A);

s = diag(S);
s = 0.5.*(s + sqrt(s.*s + 2./(rho/n)));

ret = U*diag(s)*U';
Theta = ret;
