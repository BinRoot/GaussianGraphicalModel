%% --------- TAKES A p-1xp MATRIX AND APPENDS DIAGONALS TO GET SQUARE MATRIX


% ---------- Let A be p-1xp matrix. Let d be a vector of p entries.
% ---------- These d elements are added to the diagonal positions
% ---------- of A to create a pxp square matrix.




function [X] = diag_construction(Xtilde,d)

%disp("diag_construction");
%disp(size(Xtilde)); % 49 x 50
%disp(size(d)); % 50 x 1

p = size(d,1);
%disp(p); % 50

X = zeros(p,p);

X(1,1) = d(1); 
X(2:p,1) = Xtilde(:,1);

for(j = 2:p-1)
X(1:j-1,j) = Xtilde(1:j-1,j);
X(j,j) = d(j);
X(j+1:p,j) = Xtilde(j:p-1,j);    
end;
X(1:p-1,p) = Xtilde(:,p); X(p,p) = d(p);

