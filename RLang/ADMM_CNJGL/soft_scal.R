soft_scal <- function(Y, lambda) {
  m <- dim(Y)[1]
  n <- dim(Y)[2]
  
  X <- matrix(0, m, n)

  for (t in 1:n) {
    if (norm(Y[,t],"2") > lambda) {
      X[,t] <- (norm(Y[,t],"2") - lambda)/norm(Y[,t],"2") * Y[,t]
    }
    else {
      X[,t] <- matrix(0,m,1)
    }
  }

  return(X)
}
