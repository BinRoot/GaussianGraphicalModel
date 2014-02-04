%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			TEST FILE FOR ADMM ALGORITHM FOR THE PNJGL FORMULATION                  		%
% Refer to the paper: "Node based learning of multiple Gaussian Graphical 					%
% models" Karthik Mohan, Palma London, Maryam Fazel, Daniela Witten,					%
% Su-In Lee.                                                                                    %
% http://arxiv.org/abs/1303.5145   				%
% for further details on the CNJGL formulation and the ADMM algorithm.					%	
%                                                                                               %
% CONTACT Karthik Mohan (karna@uw.edu) for any questions or comments on the code.			%
% LAST UPDATE: 8/1/2013                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



addpath ADMM_CNJGL;
addpath data_generation;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------- DEFINE PROBLEM PARAMETERS ---------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem_parameters;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- INPUT DATA --------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p = 295;
data_parameters; % Define parameters for generating simulated data

%load('./data/X.txt');
%S1 = cov(X');
%S1 = S1(1:100,1:100);
%load('./data/Y.txt');
%S2 = cov(Y');
%S2 = S2(1:100,1:100);

% Generate simulated data
[Theta_true_1, Theta_true_2, ind_m_common, ind_m_pert] = data(p);

[S1, S2] = generate_samples(Theta_true_1, Theta_true_2, n1, n2); % Replace this 
% line by user-defined sample covariance matrices S1, S2 if desired.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- RUN ADMM ALGORITHM -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic();
[Theta_1,Theta_2,V_1,V_2,iter_adm,relError] = ADM_CNJGL(S1,S2,lambda_1,lambda_2, n1, n2);
t2 = toc();

fprintf('Run time of algorithm = %f seconds \n', t2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- PLOTS AND RESULTS -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ 0. Relative error between estimate and true inverse covariance matrices

rel_err = max ( norm(Theta_true_1 - Theta_1, 'fro')/norm(Theta_true_1,'fro'), norm(Theta_true_2 - Theta_2,'fro')/norm(Theta_2, 'fro'));

fprintf('Relative error of the estimated solution = %f \n', rel_err);

% ------ 1. Compare estimated cohub network to true network

figure;
subplot(2,2,1);
imagesc(Theta_true_1 - diag(diag(Theta_true_1)));
colorbar
title('True network 1');

subplot(2,2,2);
imagesc(Theta_1 - diag(diag(Theta_1)));
title('Estimated network 1');
colorbar


subplot(2,2,3);
imagesc(Theta_true_2 - diag(diag(Theta_true_2)));
colorbar
title('True network 2');

subplot(2,2,4);
imagesc(Theta_2 - diag(diag(Theta_2)));
title('Estimated network 2');
colorbar();
