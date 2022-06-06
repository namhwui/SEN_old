library(R.utils)
library(Rfast)
library(e1071)

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
  })
  
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
    
    lapply(1:ncol(wt), function(g) {
      mu <- colSums(sweep(data, 1, wt[, g], "*")) / sum(wt[, g])
      
      s0    <- cov.wt(data, wt[, g], method = "ML")$cov
      svd_s <- svd(s0, nu = q[g], nv = 0)
      
      gamma  <- svd_s$u
      xi     <- sqrt(svd_s$d[1:q[g]])
      
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




init_hc <- function(data, G, ...) {
  hc_call <- R.utils::doCall(hclust, d = dist(data), ...)
  R.utils::doCall(cutree, tree = hc_call, k = G, ...)
}

init_kmeans <- function(data, G, ...) {
  R.utils::doCall(kmeans, x = data, centers = G, ...)$cluster
}



convergence_aitken <- function(eps = 1e-2) {
  function(model) {
    if (model$ii < 3) return(F)
    
    three <- model$logl[c(-2, -1, 0) + model$ii]
    a <- (three[3] - three[2]) / (three[2] - three[1])
    logl_asymp <- three[2] + (three[3] - three[2]) / (1 - a)
    
    if (abs(logl_asymp - three[2]) < eps) return(T)
    else return(F)
  }
}


criterion_BIC <- function(model) {
  G <- length(model$prop)
  p <- ncol(model$data)
  q <- model$q
  logl <- tail(model$logl, 1)
  k <- (G - 1) + 2*G*p + sum(p*q - q*(q+1)/2)
  2*logl - k*log(nrow(model$data))
}


EM_label <- function(model) {
  wt <- model$wt
  if (ncol(wt) == 1) return(rep(1, nrow(wt)))
  
  sapply(1:nrow(wt), function(ii) {
    which.max(wt[ii, ])
  })
}


EM_control <- function(init_label = init_kmeans, 
                       criterion = criterion_BIC, 
                       convergence = convergence_aitken(), 
                       iter = NULL, 
                       iter_max = 500) {
  list(init_label = init_label,
       criterion = criterion,
       convergence = convergence,
       iter = iter,
       iter_max = iter_max)
}


EM_until_converge <- function(model, Estep, Mstep, control) {
  
  iter_max    <- control$iter_max
  criterion   <- control$criterion
  convergence <- control$convergence
  
  model$ii   <- 1
  model$logl <- numeric(iter_max)
  
  while (model$ii <= iter_max) {
    
    model$wt  <- Estep(model)
    model     <- Mstep(model)
    model$logl[model$ii] <- Estep(model, logl.only = T)
    
    if (convergence(model)) {
      break
    } else {
      model$ii <- model$ii + 1
    }
    
    
  }
  
  model$logl      <- tail(model$logl[1:min(model$ii, iter_max)], 1)
  model$criterion <- criterion(model)
  model$label     <- EM_label(model)
  model$G <- length(model$prop)
  model$ii <- NULL
  model
}


EM_fixed_iter <- function(model, Estep, Mstep, control) {
  
  iter <- control$iter
  logl <- numeric(iter)
  
  for (ii in 1:iter) {
    model$wt <- Estep(model)
    model    <- Mstep(model)
    logl[ii] <- Estep(model, logl.only = T)
  }
  
  model$logl <- tail(logl, 1)
  model$criterion <- control$criterion(model)
  model$label     <- EM_label(model)
  model$G <- length(model$prop)
  model
}


SEN_EM <- function(data, scale = T, G = 1:6, 
               q = 1, alpha = 0.5, rho = 1,
               control = EM_control(), ...) {
  
  if (scale) data <- scale(data, center = T, scale = T)
  init_label <- control$init_label
  
  
  all_models <- lapply(G, function(g) {
    
    model <- list()
    model$data  <- data
    model$label <- init_label(data, g)
    model$prop  <- table(model$label) / length(model$label)
    model$wt    <- sapply(1:g, function(ii) model$label == ii) * 1
    model$param <- init_param(q, alpha, rho)(data, model$wt)
    
    if (is.null(control$iter)) {
      model <- EM_until_converge(model, Estep, Mstep, control)
    } else {
      model <- EM_fixed_iter(model, Estep, Mstep, control)
    }
  })
  
  all_models[[which.max(sapply(all_models, function(x) x$criterion))]]
}


ARI <- function(lab1, lab2, digits = 3, print_table = F) {
  if (length(lab1) != length(lab2)) {
    stop("Length of the label sets are not equal.")
  }
  tab <- table(lab1, lab2)
  if (print_table) {
    print(tab)
  }
  round(e1071::classAgreement(tab)$crand, digits)
}
