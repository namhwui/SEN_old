# EM steps for SEN


SEN_major <- function(param, alpha = 0.5, eps = 1e-4) {
  x <- param$gamma
  A <- 1 / (2 * sqrt(x^2 + eps))
  W <- 1 / (2 * sqrt(rowSums(x^2) + eps))
  
  SP1 <- t(A * x)
  SP2 <- t(apply(x, 1, W - max(W), "*"))
  
  alpha*SP1 + (1 - alpha)*SP2
}

FA_btS <- function(param, wt) {
  temp <- param$lambda %*% t(param$lambda)
  diag(temp) <- diag(temp) + param$psi
  beta <- t(param$lambda) %*% solve(temp)
  
  S <- cov.wt(param$data, wt = param$wt, method = "ML")$cov
  
  theta <- -beta %*% param$lambda + beta %*% S %*% t(beta)
  diag(theta) <- diag(theta) + 1
  
  list(beta = beta, S = S, theta = theta)
}

FA_midstep <- function(param, btS) {
  beta <- btS$beta
  S <- btS$S
  theta <- btS$theta
  
  A <- sweep(beta %*% S, 2, 1/param$psi, "*") * param$xi
  B <- sweep(theta, 2, param$xi, "*") * param$xi
  M <- t(sweep(param$gamma, 1, 1/param$psi, "*")) %*% param$gamma
  N <- as.matrix(diag(sweep(beta %*% S, 2, 1/param$psi) %*% param$gamma))
  
  FA_major <- A - 0.5 * t(apply(param$gamma %*% B, 1, 1/param$psi - min(param$psi), "*"))
  
  list(FA_major = FA_major, A = A, B = B, M = M, N = N, beta = beta, S = S, theta = theta)
}



Mstep_SEN <- function(model) {
  
  wt <- model$wt
  model$prop <- colSums(wt) / nrow(model$data)
  for (g in 1:model$G) {
    btS <- FA_btS(model$param[[g]], wt[, g])
    mid <- FA_midstep(model$param[[g]], btS)
    SEN <- SEN_major(param[[g]], model$alpha[g])
    
    param$mu <- colSums(sweep(model$data, 1, wt[, g], "*")) / sum(wt[, g])
    param$psi <- diag(mid$S - 2 * mid$S %*% t(mid$beta) %*% param$lambda + param$lambda %*% mid$theta %*% t(param$lambda))
    param$xi <- solve(mid$theta * mid$M) %*% mid$N
    param$gamma <- with(svd(FA_major - model$rho[g] * SEN), v %*% t(u))
    param$lambda <- sweep(param$gamma, 2, param$xi, "*")
    param$sigma <- tcrossprod(param$lambda)
    diag(param$sigma) <- diag(param$sigma) + param$psi
  }
  
  model
}

Estep_SEN <- function(model) {
  comp_dens <- sapply(1:model$G, function(g) {
    model$prop[g] * Rfast::rmvnorm(model$data, model$param[[g]]$mu, model$param[[g]]$sigma)
  })
  model$wt <- comp_dens / rowSums(comp_dens)
  model
}