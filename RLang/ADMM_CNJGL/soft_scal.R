soft_scal <- function(Y, lambda) {
  m = dim(Y)[1]
  n = dim(Y)[2]

  for (t in 1:n) {
    if (norm(Y[,t],"2") > lambda) {
      X[,t] = (norm(Y[,t],"2") - lambda)/norm(Y[,t],"2") %*% Y[,t]
    }
  }
}
