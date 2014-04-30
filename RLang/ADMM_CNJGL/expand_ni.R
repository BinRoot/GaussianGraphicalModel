expand_ni <- function(B, S, rho, n) {
#  print("expand_ni input: ")

  A <- B - S * n / (2*rho)
  e <- eigen(A)
  U <- e$vectors
  s <- e$values

#  U <- U[,ncol(U):1]
#  s <- rev(s)

#  print("U:")
#  print(U)

#  print("s:")
#  print(s)
#  readline("expanding...")

  s <- 0.5 * (s + sqrt(s*s + 2/(rho/n)))
  Theta <- U %*% diag(s) %*% t(U)


#  print("expand_ni output:")
#  print(Theta)
#  readline("press any key to continue...")
  return(Theta)
}
