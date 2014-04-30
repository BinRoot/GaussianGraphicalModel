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

% p = 295;
p = 50;
data_parameters; % Define parameters for generating simulated data


load('./data/Y.txt');
Y = Y(:,2:2);


% X = importdata('../data/X.txt');
load('../data/X.txt');
#save -6 X.mat X
disp(size(X))


% all 1 value entries
X1 = X(find(Y+1), :);

% all -1 value entries
X2 = X(find(Y-1), :);


[x1rows, x1cols] = size(X1);
[x2rows, x2cols] = size(X2);

[sorted, indices] = sort(std(X));
X1 = X1(1:x1rows, indices(end-9:end));
X2 = X2(1:x2rows, indices(end-9:end));

% Input data
%X = [0.05, -0.04, 0.06, 0.02, -0.01;
%     0.02, 0.01, -0.06, 0.01, -0.04;
%     0.05, 0.04, -0.02, 0.01, -0.01]

S1 = cov(X1);

disp("cov matrix is of size")
disp(size(S1))

p=size(S1)(1);

%load('./data/Y.txt');
#save -6 Y.mat Y
%Y = Y(:,2:2);
%disp(size(Y));
%S2 = cov(Y);
%disp(size(S2));
%S2 = S2(1:10,1:10);

S2 = cov(X2);

#disp("input is")
#disp(S1)

% Generate simulated data
%[Theta_true_1, Theta_true_2, ind_m_common, ind_m_pert] = data(p);
%[S1, S2] = generate_samples(Theta_true_1, Theta_true_2, n1, n2); % Replace this 
% line by user-defined sample covariance matrices S1, S2 if desired.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- RUN ADMM ALGORITHM -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic();
[Theta_1,Theta_2,V_1,V_2,iter_adm,relError] = ADM_CNJGL(S1,S1,lambda_1,lambda_2, n1, n2);
t2 = toc();

fprintf('Run time of algorithm = %f seconds \n', t2);

disp(Theta_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- PLOTS AND RESULTS -------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ 0. Relative error between estimate and true inverse covariance matrices

% rel_err = max ( norm(Theta_true_1 - Theta_1, 'fro')/norm(Theta_true_1,'fro'), norm(Theta_true_2 - Theta_2,'fro')/norm(Theta_2, 'fro'));

% fprintf('Relative error of the estimated solution = %f \n', rel_err);

% ------ 1. Compare estimated cohub network to true network

%figure;
%subplot(2,2,1);
%imagesc(Theta_true_1 - diag(diag(Theta_true_1)));
%colorbar
%title('True network 1');

%subplot(2,1,1);
%imagesc(Theta_1 - diag(diag(Theta_1)));
%title('Estimated network 1');
%colorbar();

res1 = Theta_1 - diag(diag(Theta_1));
res1 = reshape (res1,p*p,1);
[res1,res1_indices] = sort(res1);
res1 = res1(end-10: end);
res1_indices = res1_indices(end-10: end);
res1_indices = [ idivide(res1_indices-1,p)+1  res1_indices-idivide(res1_indices-1,p)*p ];


save res1.mat res1;
save res1_indices.mat res1_indices;
disp("data X values saved");

%subplot(2,2,3);
%imagesc(Theta_true_2 - diag(diag(Theta_true_2)));
%colorbar
%title('True network 2');

%subplot(2,1,2);
%imagesc(Theta_2 - diag(diag(Theta_2)));
%title('Estimated network 2');
%colorbar();

res2 = Theta_2 - diag(diag(Theta_2));
res2 = reshape (res2,p*p,1);
[res2,res2_indices] = sort(res2);
res2 = res2(end-10: end);
res2_indices = res2_indices(end-10: end);
res2_indices = [ idivide(res2_indices-1,p)+1  res2_indices-idivide(res2_indices-1,p)*p ];

disp("data Y values saved");
save res2.mat res2;
save res2_indices.mat res2_indices;
