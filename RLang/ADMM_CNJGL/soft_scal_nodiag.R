soft_scal_nodiag <- function(Y,lambda) {
  p <- dim(Y)[2];
  X <- matrix(0,p,p);
  Xtilde = matrix(0,p-1,p-1);
  Ytilde = nodiag_construction(Y);
  Xtilde = soft_scal(Ytilde,lambda);
  X = diag_construction(Xtilde,diag(Y));
}
