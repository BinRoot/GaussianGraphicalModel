%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------- ADMM algorithm for the PNJGL formulation ------------ %%

% Refer to the paper: "Structured learning of Gaussian Graphical 
% models" Karthik Mohan, Mike Chung, Seungeyop Han, Daniela Witten,
% Su-In Lee and Maryam Fazel.
% http://students.washington.edu/karna/Papers/nips12_camera_ready.pdf
% for details.

% -------- CONTACT Karthik Mohan (karna@uw.edu) for any questions or
% -------- comments on the code -------------------------------------%
% ------------------ LAST UPDATE: 2/28/2013 --------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Theta_1,Theta_2,V] = ADM_PNJGL(S1,S2,lambda_1,lambda_2, n_1, n_2) 
% Initialize

n = size(S1,1); % number of variables  

%------ Initialization -------------- %

 % Primal variables
eye_n = eye(n);
Theta_1 = eye_n; 
Theta_2 = eye_n; 
Theta_1_old = eye_n; 
Theta_2_old = eye_n;
Z1 = eye_n; 
Z2  = eye_n;

 % Dual Variables
V = zeros(n,n); W = zeros(n,n); Lambda = zeros(n,n); Gamma = zeros(n,n);
Q1 = zeros(n,n); Q2 = zeros(n,n); 


% ----------- Set parameters for ADMM ---------- %
ADM_parameters_PNJGL; 



%------------------------------
% BEGIN ALGORITHM
%------------------------------


rho = opts.rho_init;


%----------------------------------------
% Begin outer loop for homotopy parameter
%----------------------------------------

for(r = 1:homotopy_size)
    
  
    %---------------------
    % Begin ADM 
    %---------------------
    
 for(i = 1:opts.maxiter)
     
     % Update 6 primal variables and 4 dual variables in the ADM formulation
     
     Theta_1 = expand_n(1/(2*rho)*(rho*(Theta_2 + V + W+ Z1) - (Q1 + (S1*n_1) + Lambda)), rho, n_1);
     Theta_2 = expand_n(1/(2*rho)*(rho*(Theta_1 - (V + W) + Z2) - (Q2 + (S2*n_2) - Lambda)), rho, n_2);
     Z1 = soft_thresh(Theta_1 + Q1/rho, lambda_1/rho);
     Z2 = soft_thresh(Theta_2 + Q2/rho, lambda_1/rho);
     V = soft_scal(0.5*(Theta_1 + W' - Theta_2 - W) + 1/(2*rho)*(Lambda - Gamma), lambda_2/(2*rho));
     W = 0.5*(Theta_1 - Theta_2 + V' - V) + 1/(2*rho)*(Lambda + Gamma');
     
     Lambda = Lambda + rho*(Theta_1 - Theta_2 - (V + W));
     Gamma = Gamma + rho*(V - W');
     Q1 = Q1 + rho*(Theta_1 - Z1);
     Q2 = Q2 + rho*(Theta_2 - Z2);
    
     
    relerr_temp = max(norm(Theta_1 - Theta_1_old,'fro')/norm(Theta_1_old,'fro'), norm(Theta_2 - Theta_2_old,'fro')/norm(Theta_2_old,'fro'));

    if(relerr_temp <= opts.eps)
        break;
    end;
    
    
    Theta_1_old = Theta_1; 
    Theta_2_old = Theta_2;
 
 end;
     %---------------------
     % END ADM 
     %---------------------
 
rho = rho * opts.homotopy;
end;


%-----------------------------
% END ALGORITHM
%-----------------------------

