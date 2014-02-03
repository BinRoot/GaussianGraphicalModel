% ------------------- GENERATE DATA (SAMPLE COVARIANCE MATRICES) FOR THE ADM ALGORITHM -------- %%

% Palma London:  update 8/20/12 
% Karthik Mohan: update 11/5/12

function [S_1, S_2] = generate_samples(Theta_true_1, Theta_true_2, n_1, n_2) 
    
    p = size(Theta_true_1,1); % number of variables 
    
    X_1 = mgd(n_1, p, zeros(p, 1), inv(Theta_true_1) ); % Generates ssize samples from True Covariance Matrix i
    S_1 = cov(X_1); % Sample Covariance matrix i

    X_2 = mgd(n_2, p, zeros(p, 1), inv(Theta_true_2) ); % Generates ssize samples from True Covariance Matrix i
    S_2 = cov(X_2); % Sample Covariance matrix i


end
