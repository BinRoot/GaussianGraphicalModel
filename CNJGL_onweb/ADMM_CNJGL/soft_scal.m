%% ------------------- SOFT-SCALING OPERATOR ----------------------- %%

% -------------------- LAST UPDATE: 4/24/2012 ------------------------%

% ---------------- minimize_X \frac{1}{2}||X - Y||_2^2 + \lambda*sum_i||X_i||_2

% ---------------- solution is  X_i = max(||Y_i||_2 - \lambda, 0)*Y_i/||Y_i||_2
% %


% -------------------- LAST UPDATE: 4/19/2012 ------------------------ %


function [X] = soft_scal(Y,lambda)

m = size(Y,1);
n = size(Y,2);


for(t = 1:n)
if(norm(Y(:,t),2) > lambda)
    X(:,t) = (norm(Y(:,t),2) - lambda)/norm(Y(:,t),2)*Y(:,t);
else
    X(:,t) = zeros(m,1);
end;
end;



