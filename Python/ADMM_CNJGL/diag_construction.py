import numpy

def diag_construction(xtilde, d):
    p = d.shape[0]
    X = numpy.zeros((p, p))
    X[0][0] = d[0]
    X[1:p][0] = xtilde[:][1]
    for j in range (1, p-1):
        X[1:j-1][j] = xtilde[1:j-1][j]
        X[j][j] = d[j]
        X[j+1:p-1][j] = xtilde[j:p-2][j]
    X[0:p-2][p-1] = xtilde[:][p-1]
    X[p-1][p-1] = d[p-1]

'''

X(1,1) = d(1); 
X(2:p,1) = Xtilde(:,1);

for(j = 2:p-1)
X(1:j-1,j) = Xtilde(1:j-1,j);
X(j,j) = d(j);
X(j+1:p,j) = Xtilde(j:p-1,j);    
end;
X(1:p-1,p) = Xtilde(:,p); 
X(p,p) = d(p);


'''

