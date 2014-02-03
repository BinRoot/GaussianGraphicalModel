%% ------------------- SOFT-SCALING OPERATOR ----------------------- %%

% -------------------- LAST UPDATE: 2/28/2013 ------------------------%

% ---------------- minimize_X \frac{1}{2}||X - Y||_2^2 + \lambda*sum_i||X_i||_2

% ---------------- solution is  X_i = max(||Y_i||_2 - \lambda, 0)*Y_i/||Y_i||_2
% %





function [X] = soft_scal(Y,lambda)

if(nargin == 1)
    disp('Not enough inputs');
    disp('Enter both Y and lambda');
    return;
end;

m = size(Y,1);
n = size(Y,2);

% Preallocate
X = zeros(m,n);


for(t = 1:n)
    
if(norm(Y(:,t),2) > lambda)
    X(:,t) = (norm(Y(:,t),2) - lambda)/norm(Y(:,t),2)*Y(:,t);
end;

end;



