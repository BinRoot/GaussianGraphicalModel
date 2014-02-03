%% ---------- REMOVE DIAGONAL ELEMENTS OF A MATRIX AND REARRANGE -------- %%

% ----------- A pxp matrix X is taken as input, its diagonals are removed.
% ----------- The resulting matrix is of size p-1xp.

% ----------- LAST UPDATE: 8/15/2012 ------------------------------- %

function Xtilde = nodiag_construction(X);

p = size(X,1);
Xtilde = zeros(p-1,p-1);


Xtilde(:,1) = X(2:p,1);
for(j = 2:p-1)
    Xtilde(:,j) = X([1:j-1,j+1:p],j);
end;
Xtilde(:,p) = X(1:p-1,p);