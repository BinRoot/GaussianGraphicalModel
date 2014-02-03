%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- The Expand Operator ----------------------------- %
% --- Expand_n(A, rho, n) = argmin_{Theta} -n\log\det(Theta) + rho*||Theta -   %
% A||_F^2								       %	
% --- Let A = U*Sigma*U' and rho1 = rho/n;				       %	
% Set Sigma2_{ii} = 0.5*(Sigma_{ii} + sqrt(Sigma_{ii}^2 + 2/rho1));            %  			 
% Expand_n(A, rho,n) = U*Sigma2*U'					       %	
%									       %	
% LAST UPDATE: 3/1/2013							       %	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Theta] = expand_n(A, rho, n)

[U,S] = eig(A); 
s = diag(S);
s = 0.5.*(s + sqrt(s.*s + 2./(rho/n)));

% a way to get around using the Identity matrix. 
Theta = U*diag(s)*U';

end
