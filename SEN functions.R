# EM steps for SEN


SEN_major <- function(param, eps = 1e-4) {
  x <- param$gamma
  A <- 1 / (2 * sqrt(x^2 + eps))
  W <- 1 / (2 * sqrt(rowSums(x^2) + eps))
  
  SP1 <- t(A * x)
  SP2 <- t(sweep(x, 1, W - max(W), "*"))
  
  param$alpha*SP1 + (1 - param$alpha)*SP2
}

SEN_update <- function(data, param, wt) {
  S      <- cov.wt(data, wt = wt, method = "ML")$cov
  invpsi <- 1/param$psi
  beta   <- with(param, {t(lambda) %*% solve(sigma)})
  
  theta       <- -beta %*% param$lambda + beta %*% S %*% t(beta)
  diag(theta) <- diag(theta) + 1
  
  M <- sweep(t(param$gamma), 2, invpsi, "*") %*% param$gamma
  N <- diag(beta %*% sweep(S, 2, invpsi, "*") %*% param$gamma)
  param$xi <- c(solve(theta * M) %*% N)
  
  A     <- sweep(beta, 1, param$xi, "*") %*% sweep(S, 2, invpsi, "*")
  B     <- sweep(sweep(theta, 2, param$xi, "*"), 1, param$xi, "*")
  mat   <- A - B %*% t(param$gamma) %*% diag(invpsi - max(invpsi))
  major <- SEN_major(param)
  
  param$gamma  <- with(svd(mat - param$rho * major), v %*% t(u))
  param$lambda <- sweep(param$gamma, 2, param$xi, "*")
  
  param$psi         <- diag(S - 2*S%*%t(param$lambda%*%beta) + param$lambda%*%theta%*%t(param$lambda))
  param$sigma       <- tcrossprod(param$lambda)
  diag(param$sigma) <- diag(param$sigma) + param$psi
  
  param$mu <- colSums(sweep(data, 1, wt, "*")) / sum(wt)
  
  param
}


Mstep <- function(model) {
  
  wt <- model$wt
  model$prop <- colSums(wt) / nrow(model$data)
  for (g in 1:ncol(wt)) {
    model$param[[g]] <- SEN_update(data, model$param[[g]], wt[, g])
  }
  
  model
}



Estep <- function(model, logl.only = F) {
  comp_dens <- sapply(1:ncol(model$wt), function(g) {
    model$prop[g] * Rfast::dmvnorm(model$data, model$param[[g]]$mu, model$param[[g]]$sigma)
    #print(cbind(model$param[[g]]$mu, colMeans(model$data[model$label == g, ])))
  })
  #print(comp_dens)
  if (logl.only) {
    return(sum(log(rowSums(comp_dens))))
  } else {
    comp_dens / rowSums(comp_dens)
  }
}

init_param <- function(q, alpha, rho) {
  function(data, wt) {
    G <- ncol(wt)
    q <- rep(q, G)
    alpha <- rep(alpha, G)
    rho <- rep(rho, G)
    #print(q)
    
    lapply(1:ncol(wt), function(g) {
      mu <- colSums(sweep(data, 1, wt[, g], "*")) / sum(wt[, g])
      
      s0    <- cov.wt(data, wt[, g], method = "ML")$cov
      svd_s <- svd(s0, nu = q[g], nv = 0)
      #print(svd_s)
      gamma  <- svd_s$u
      xi     <- sqrt(svd_s$d[1:q[g]])
      #print(xi)
      lambda <- sweep(gamma, 2, xi, "*")
      psi    <- abs(diag(s0 - lambda %*% t(lambda)))
      
      sigma       <- lambda %*% t(lambda)
      diag(sigma) <- diag(sigma) + psi
      
      list(mu     = mu,
           gamma  = gamma,
           xi     = xi,
           lambda = lambda,
           psi    = psi,
           sigma  = sigma,
           alpha  = alpha[g],
           rho    = rho[g])
    })
    
  }
}
