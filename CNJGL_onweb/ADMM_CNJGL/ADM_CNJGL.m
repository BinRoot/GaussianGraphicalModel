%%--------------- ADMM algorithm for CNJGL -------------- %%
% ------------------- KARTHIK MOHAN, UW ------------------ %

% ------------------------ LAST UPDATE: 8/1/2013 ---------------------- %



function [Theta_1,Theta_2,V_1,V_2,iter_adm, relError] = ADM_CNJGL(S_1,S_2,lambda_1,lambda_2,n_1,n_2)


%------------------------------
% BEGIN ALGORITHM
%------------------------------
p = size(S_1,1);

%----------------------
% INITIALIZE VARIABLES
%----------------------
Theta_1_old = eye(p); Theta_2_old = eye(p);
Theta_1 = eye(p); Theta_2 = eye(p);
V_1 = Theta_1_old/2; V_2 = Theta_2_old/2;
W_1 = V_1'; W_2 = V_2'; 
Z_1 = Theta_1; Z_2 = Theta_2;
Gamma_1 = zeros(p,p); Gamma_2 = Gamma_1; 
Lambda_1 = zeros(p,p); Lambda_2 = zeros(p,p); 
Q_1 = zeros(p,p); Q_2 = Q_1; 
%}
%----------------------
% LOAD PARAMETERS 
%----------------------
ADM_parameters;
rho = opts.rho_init;
iter_adm = zeros(opts.homotopy_size,1);

relError = zeros(opts.maxiter, 1);

%----------------------------------------
% Begin outer loop for homotopy parameter
%----------------------------------------

disp ("homotopy_size")
disp (opts.homotopy_size)

disp ("maxiter")
disp (opts.maxiter)

for (r = 1:opts.homotopy_size)
    
    %---------------------
    % Begin ADM 
    %---------------------
 for (i = 1:opts.maxiter)
     %Update 8 primal variables and 6 dual variables in the ADM formulation

     % Update of Theta_1 and Theta_2
     Theta_1 = expand_ni(1/(2*rho)*(rho*(V_1 + W_1 + Z_1) - (Gamma_1 + Lambda_1)),S_1,rho,n_1);
     Theta_2 = expand_ni(1/(2*rho)*(rho*(V_2 + W_2 + Z_2) - (Gamma_2 + Lambda_2)),S_2,rho,n_2);

     % Update of Z_1 and Z_2
     Z_1 = soft_thresh(Theta_1 + Lambda_1/rho, lambda_1/rho);
     Z_2 = soft_thresh(Theta_2 + Lambda_2/rho, lambda_1/rho);

     % UPDATE OF V_1 and V_2
     C_1 = 1/(2*rho)*(rho*(W_1' + Theta_1 - W_1) - (Q_1 - Gamma_1));
     C_2 = 1/(2*rho)*(rho*(W_2' + Theta_2 - W_2) - (Q_2 - Gamma_2));

%     disp("C_1")
%     disp(size(C_1))
     d_1 = diag(C_1); d_2 = diag(C_2);
     N_1 = nodiag_construction(C_1);
     N_2 = nodiag_construction(C_2);
%     disp("N_1")
%     disp(size(N_1))
     N = [N_1; N_2];
%     disp("N")
%     disp(size(N))
     H = soft_scal(N,lambda_2/(2*rho)); d = [d_1;d_2];

%     disp("d_1")
%     disp(size(d_1))
%     disp("H")
%     disp(size(H))

%     disp("diagcon H")
%     disp(size(H(1:p-1,:)))

     V_1 = diag_construction(H(1:p-1,:),d_1);

%     disp("diagcon2 H")
%     disp(size(H(p:2*p-2,:)))

     V_2 = diag_construction(H(p:2*p-2,:),d_2);

     % Update of W_1 and W_2
     W_1 = 1/(2*rho)*(rho*(V_1' + Theta_1 - V_1) + Gamma_1 + Q_1');
     W_2 = 1/(2*rho)*(rho*(V_2' + Theta_2 - V_2) + Gamma_2 + Q_2');

     % Update of Dual variables
     Gamma_1 = Gamma_1 + rho*(Theta_1 - V_1 - W_1);
     Gamma_2 = Gamma_2 + rho*(Theta_2 - V_2 - W_2);
     Lambda_1 = Lambda_1 + rho*(Theta_1 - Z_1);
     Lambda_2 = Lambda_2 + rho*(Theta_2 - Z_2);
     Q_1 = Q_1 + rho*(V_1 - W_1');
     Q_2 = Q_2 + rho*(V_2 - W_2');

%     disp("i")
%     disp(i)
%     disp("r")
%     disp(r)
%     disp("relError index")
%     disp(i + ((r - 1)*opts.maxiter))

     relError(i + ((r - 1)*opts.maxiter), 1) = ...
         max(norm(Theta_1 - Theta_1_old,'fro')/norm(Theta_1,'fro'), ...
             norm(Theta_2 - Theta_2_old,'fro')/norm(Theta_2,'fro'));

    if(relError(i + ((r - 1)*opts.maxiter), 1) <= opts.eps)
        break;
    end;
    Theta_1_old = Theta_1; Theta_2_old = Theta_2;
 end;
     %---------------------
     % END ADM 
     %---------------------
 
rho = rho*opts.homotopy;
iter_adm(r,1) = i;
end;


%-----------------------------
% END ALGORITHM
%-----------------------------
