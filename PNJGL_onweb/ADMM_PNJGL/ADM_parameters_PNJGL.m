%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------- SET PARAMETERS FOR ADM ALGORITHM -------------- %
%								%
% 		LAST UPDATE: 3/1/2013				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




opts.maxiter = 1000;
opts.eps = 1e-4; % Tolerance for convergence of the ADMM algorithm
opts.rho_init = 1;
opts.homotopy = 5; % Scaling parameter for rho
opts.rho_max = 5;

if(opts.rho_max/opts.rho_init == 1)
    homotopy_size = 1;
else
    homotopy_size = log(opts.rho_max)/log(opts.homotopy) + 1;
    homotopy_size = floor(homotopy_size);
end;
