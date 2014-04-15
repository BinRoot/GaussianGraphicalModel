## ---------- REMOVE DIAGONAL ELEMENTS OF A MATRIX AND REARRANGE -------- %%

# ----------- A pxp matrix X is taken as input, its diagonals are removed.
# ----------- The resulting matrix is of size p-1xp.

# ----------- LAST UPDATE: 8/15/2012 ------------------------------- %

nodiag_construction <- function(X) {
  p <- dim(X)[1]
  Xtilde <- matrix(0,p-1,p)
  Xtilde[,1] <- X[2:p,1]
  for (j in 2:(p-1)) {
    Xtilde[,j] <- c(X[1:(j-1),j], X[(j+1):p,j])
  }
  Xtilde[,p] <- X[1:(p-1),p];
  return(Xtilde)
}

