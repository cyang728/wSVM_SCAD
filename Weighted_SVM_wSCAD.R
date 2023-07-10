scadsvc <- function(lambda1 = 0.01, x, y, a = 3.7, tol = 10^(-4), sample.weights = w_hat,
                    seed = 123, maxIter = 1, verbose = FALSE){
  # SCAD weigghted SVM classification
  #
  # Input:
  #   xtrain : n-by-d data matrix to train (n chips/patients, d clones/genes, d>>n )
  #   ytrain : column vector of target {-1, 1}'s (for n chips/patiens )
  #   a : tuning parameter in scad function (default: 3.7 or whatever the paper uses)
  #   lambda1 : tuning parameter in scad function (default : 2 or whatever the paper uses)
  #   tol: the cut-off value to be taken as 0
  #
  # Cheng-Han Yang 
  # The R package penalizedSVM, developed by Axel Benner (2007/03/15), is used to modify the code.
  
  # SVM mit variable selection (clone selection) using scad penalty.
  
  # checks
  if (nrow(x) != length(y)) stop("Wrong dimensions: nrow(x) should be equal to length(y) !")
  if (nlevels(as.factor(y)) != 2) stop(paste("We need 2 classes, currently have:", paste(levels(as.factor(y)), collapse=", ")))
  
  xtrain <- x
  # change class labels to -1 and 1, y traqin is a numeric vector !
  ytrain <- 2 * as.numeric(as.factor(y)) - 3
  
  # start with linear svm:
  require(e1071, quietly = TRUE)
  require(WeightSVM, quietly = TRUE) # for weighted SVM
  require(corpcor, quietly = TRUE) # for pseudoinverse
  require(statmod, quietly = TRUE) # for matrix multplications vecmat and matvec
  
  set.seed(seed)
  dat_tmp <- cbind(xtrain, ytrain)
  colnames(dat_tmp)[ncol(xtrain) + 1] <- "trt"
  obj <- WeightSVM::wsvm(trt ~ ., weight = sample.weights, type = "C-classification",
                         data = dat_tmp, kernel = "linear")
  index <- obj$index
  # type of svm
  type <- obj$type
  # w: coefficients times the training labels * SV
  w <- apply(as.vector(obj$coefs) * obj$SV, 2, sum)
  
  # rho = Bias: A scalar value representing the bias or threshold of the SVM
  # classifier,  which is the negative intercept.
  b <- obj$rho
  diff <- 1000 # was 1, in original script also 1000 (last modification: 18-March-2009)
  ntrain <- sum(sample.weights)
  d <- ncol(xtrain)
  xind <- 1:d
  i <- 1
  if (verbose) print("start iterations:")
  while (diff > tol) {
    # should write i in the same position
    # cat("\b\b\b\b\b",i)
    # flush.console()
    
    if (!is.null(maxIter)) {
      if (i > maxIter) {
        if (verbose) print(paste("max. iteration number of", maxIter, "has been reached. Stop iterations "))
        nw <- w
        nb <- b
        break
      }
    }
    
    x1 <- cbind(rep(1, nrow(xtrain)), xtrain)
    x1 <- as.matrix(x1)
    sgnres <- as.vector(ytrain - x1 %*% c(b, w))
    # important!!!! : sometimes a point is lying exactly on the hyperline --> sgnres = 0
    # produce errors --> move this randomly at the one or the other size.
    sgnres[sgnres == 0] <- sample(c(1, -1), 1) * 10^(-100)
    res <- abs(sgnres)
    y0 <- ytrain / res
    
    # ###
    # D = 1/(2*ntrain) * diag(1/res)
    # ###
    # save as a vector D_vec
    D_vec <- wt / (2 * ntrain) * (1 / res)
    aw <- abs(w)
    dp <- lambda1 * (aw <= lambda1) + (a * lambda1 - aw) / (a - 1) * (lambda1 < aw & aw <= a * lambda1)
    Q1_vec <- c(0, dp / aw)
    
    P <- 0.5 * t(ytrain + y0) %*% diag(wt) %*% x1 / ntrain
    
    # ###
    # Q1 = diag( c(0, dp/aw))
    # Q = t(x1) %*% D %*% x1 + Q1
    # ps_Q <- pseudoinverse(Q)
    # ###
    
    # inv_Q is sometimes too large
    # inv_Q <- .find.inverse(U = t(x1), D_vec = D_vec, A_vec = Q1_vec)
    # nwb = inv_Q %*% t(P)
    
    # nwb = inv( t(x1) %*% D %*% x1 + Q1) %*% t(P)
    # -->
    nbw <- .calc.mult.inv_Q_mat3(U = t(x1), D_vec = D_vec, A_vec = Q1_vec, mat2 = t(P), n.thr = 10^4)
    # summary(nbw)
    
    nw <- nbw[-1]
    nb <- nbw[1]
    diff <- sqrt(sum((nbw - c(b, w))^2))
    ind <- abs(nw) > 0.001 # oracle threshold
    if (sum(!is.na(ind)) > 0 & sum(ind) > 0) {
      w <- nw[ind]
      xtrain <- xtrain[, ind, drop = FALSE]
      xind <- xind[ind]
      b <- nb
    } else {
      diff <- tol / 2
      xind <- 0
    }
    i <- i + 1
  }
  if (verbose) print(paste("scad converged in", i - 1, "iterations"))
  
  ind <- abs(nw) > 0.001 # oracle threshold
  
  if (sum(!is.na(ind)) > 0 & sum(ind) > 0) {
    
    w <- nw[ind]
    names(w) <- colnames(xtrain)[ind]
    b <- nb
    f <- as.vector(as.matrix(xtrain) %*% w + b)
    
    # xqx = 0.5 * x1 %*% inv_Q %*% t(x1)  - we don't have a huge quadratic matrix inv_Q, use the same trick as for nwb
    qx <- .calc.mult.inv_Q_mat3(U = t(x1), D_vec = D_vec, A_vec = Q1_vec, mat2 = t(x1), n.thr = 10^4)
    xqx <- 0.5 * as.matrix(x1) %*% qx
    
    BIC <- sum(wt * pmax(0, 1 - ytrain * f)) + length(xind) * log(sum(wt))
    
    # Output:
    #   w : direction vector of separating hyperplane
    #   b : the 'bias'
    #   xind : Indices of remained variables
    #   fitted : Fit of regularized SVM (for all patients with reduced set of genes )
    ret <- list(w = w, b = b, xind = xind, index = index, xqx = xqx, fitted = f,
                type = type, lambda1 = lambda1, iter = i - 1, BIC = BIC)
    
    class(ret) <- c("scadsvm", "penSVM")
    return(ret)
  } else {
    return(-1)
  }
}

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
