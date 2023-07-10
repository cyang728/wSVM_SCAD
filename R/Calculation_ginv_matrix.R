.calc.mult.inv_Q_mat3 <- function(U, D_vec, A_vec, mat2 = NULL, n.thr = 500) {
  if (is.null(mat2)) mat2 <- diag(1, nrow = nrow(U), ncol = nrow(U))
  
  if (nrow(U) < n.thr) {
    Q <- U %*% diag(D_vec) %*% t(U) + diag(A_vec)
    
    if (rank.condition(Q)$rank == nrow(Q) & nrow(Q) == ncol(Q)) {
      inv_Q <- solve(Q)
    } else {
      inv_Q <- pseudoinverse(Q)
    }
    
    if (all(is.na(inv_Q))) {
      nbw <- .calc.mult.inv_Q_mat2(U = U, D_vec = D_vec, A_vec = A_vec, mat2 = mat2, n.thr = 0)
    } else {
      nbw <- inv_Q %*% mat2
    }
  } else {
    if (all(A_vec == 0)) {
      A <- diag(A_vec)
      D <- diag(D_vec)
      inv_U <- pseudoinverse(U)
      
      nbw <- t(inv_U) %*% (ginv(D_vec) %*% (inv_U %*% mat2))
    } else {
      eps <- 10^-8
      A1_vec <- A_vec
      A1_vec[A1_vec == 0] <- eps
      inv_A1_vec <- 1 / A1_vec
      
      nbw.term1 <- vecmat(inv_A1_vec, mat2)
      
      tmp.mat <- ginv(diag(D_vec)) + matvec(t(U), inv_A1_vec) %*% U
      
      if (rank.condition(tmp.mat)$rank == nrow(tmp.mat) & nrow(tmp.mat) == ncol(tmp.mat)) {
        try(tt1 <- solve(tmp.mat))
      }
      
      if (!exists("tt1")) try(tt1 <- pseudoinverse(tmp.mat))
      
      if (!exists("tt1")) try(tt1 <- ginv(tmp.mat))
      
      if (!exists("tt1")) stop("Error: can not calculate inverse matrix for 'diag(1/D_vec) + matvec(t(U), inv_A1_vec )%*% U'")
      
      nbw.term2 <- vecmat(inv_A1_vec, U) %*% (tt1 %*% (matvec(t(U), inv_A1_vec) %*% mat2))
      
      nbw <- nbw.term1 - nbw.term2
    }
  }
  
  return(nbw)
}
