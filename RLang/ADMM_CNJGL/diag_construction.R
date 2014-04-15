## --------- TAKES A p-1xp MATRIX AND APPENDS DIAGONALS TO GET SQUARE MATRIX
# ---------- Let A be p-1xp matrix. Let d be a vector of p entries.
# ---------- These d elements are added to the diagonal positions
# ---------- of A to create a pxp square matrix.

diag_construction <- function(xtilde, d) {
  p <- length(d)
  x <- matrix(0, p, p)
  x[1,1] <- d[1]
  x[2:p,1] <- xtilde[,1]

  for (j in 2:p-1) {
    print(j)
    x[1:(j-1),j] <- xtilde[1:(j-1),j]
    x[j,j] <- d[j]
    print("--")    
    print((j+1):p)
    print(x[(j+1):p,j])
    print(xtilde[j:(p-1),j])
    print("--")
    x[(j+1):p,j] <- xtilde[j:(p-1),j]
  }

  x[1:(p-1),p] <- xtilde[,p]
  x[p,p] <- d[p]

  return(x)                                
}
