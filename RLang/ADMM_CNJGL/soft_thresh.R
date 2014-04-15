soft_thresh <- function(Y, lambda) {
  b <- abs(Y) - lambda
  for (c in 1:dim(b)[2]) {
    for (r in 1:dim(b)[1]) {
      if (b[r,c] < 0) {b[r,c] = 0}                  
    }
  }
  return(sign(Y) * b)
}
