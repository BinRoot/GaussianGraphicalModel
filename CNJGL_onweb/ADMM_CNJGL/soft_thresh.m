%% -------------------- SOFT-THRESHOLDING OPERATOR ------------------- %%


% ---------------- minimize 1/2*||X - Y||_2^2 + lambda*||X||_1 -------- %



% ----------------- LAST UPDATE: 4/19/2012 ---------------------------- %



function [X] = soft_thresh(Y,lambda)

disp("soft_thresh:")
disp(size(Y))

X = Y;
if(nargin == 1)
    disp('Not enough inputs');
    disp('Enter both Y and lambda');
    return;
    %lambda = 1;
end;
        
X = sign(Y).*(max(abs(Y) - lambda,0));
