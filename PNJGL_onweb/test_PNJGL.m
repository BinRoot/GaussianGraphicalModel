%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			TEST FILE FOR ADMM ALGORITHM FOR THE PNJGL FORMULATION                      %
%  Refer to the paper: "Node based learning of multiple Gaussian Graphical          			%
% models" Karthik Mohan, Palma London, Maryam Fazel, Daniela Witten, 
% Su-In Lee.                                                                                     %
% http://arxiv.org/abs/1303.5145                                                        %
% for further details on the PNJGL formulation and the ADMM algorithm.                      %	
%                                                                                       %
% CONTACT Karthik Mohan (karna@uw.edu) for any questions or comments on the code.			%
% LAST UPDATE: 8/1/2013                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



addpath ADMM_PNJGL;
addpath data_generation;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------- DEFINE PROBLEM PARAMETERS ---------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem_parameters;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- INPUT DATA --------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_parameters; % Define parameters for generating simulated data
p = 295;
% Generate simulated data
load('./data/X.txt');
S1 = cov(X');
%S1 = S1(1:100,1:100);
load('./data/Y.txt');
S2 = cov(Y');

[Theta_true_1, Theta_true_2, ind_m_pert] = data(p, ...
umin_sparse,umax_sparse, umin_pert, umax_pert, sparsity_prob, m_pert);

% [S1, S2] = generate_samples(Theta_true_1, Theta_true_2, n1, n2); % Replace this 
% line by user-defined sample covariance matrices S1, S2 if desired.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- RUN ADMM ALGORITHM -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = tic;
[Theta_1,Theta_2,V] = ADM_PNJGL(S1,S2,lambda_1,lambda_2, n1, n2);
t2 = toc(t1);

fprintf('Run time of algorithm = %f seconds \n', t2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- PLOTS AND RESULTS -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ 0. Relative error between estimate and true inverse covariance matrices

rel_err = max ( norm(Theta_true_1 - Theta_1, 'fro')/norm(Theta_true_1,'fro'), norm(Theta_true_2 - Theta_2,'fro')/norm(Theta_2, 'fro'));

fprintf('Relative error of the estimated solution = %f \n', rel_err);

% ------ 1. Compare estimated networks to true networks

figure;
subplot(2,2,1);
imagesc(Theta_true_1 - diag(diag(Theta_true_1)));
colorbar
title('True network 1');

subplot(2,2,2);
imagesc(Theta_1 - diag(diag(Theta_1)));
title('Estimated network 1');
colorbar

res1 = Theta_1 - diag(diag(Theta_1));
res1 = reshape (res1,87025,1);
[res1,res1_indices] = sort(res1);
res1 = res1(end-10: end);
res1_indices = res1_indices(end-10: end);
res1_indices = [ idivide(res1_indices-1,295)+1  res1_indices-idivide(res1_indices-1,295)*295 ];

disp("data X values:");
disp(res1);
disp(res1_indices);


subplot(2,2,3);
imagesc(Theta_true_2 - diag(diag(Theta_true_2)));
colorbar
title('True network 2');

subplot(2,2,4);
imagesc(Theta_2 - diag(diag(Theta_2)));
title('Estimated network 2');
colorbar

res2 = Theta_2 - diag(diag(Theta_2));
res2 = reshape (res2,87025,1);
[res2,res2_indices] = sort(res2);
res2 = res2(end-10: end);
res2_indices = res2_indices(end-10: end);
res2_indices = [ idivide(res2_indices-1,295)+1  res2_indices-idivide(res2_indices-1,295)*295 ];

disp("data Y values:");
disp(res2);
disp(res2_indices);

