library(R.utils)

hc_label <- function(data, G, ...) {
  hc_call <- R.utils::doCall(hclust, d = dist(data), ...)
  R.utils::doCall(cutree, tree = hc_call, k = G, ...)
}

kmeans_label <- function(data, G, ...) {
  R.utils::doCall(kmeans, x = data, centers = G, ...)$cluster
}



convergence_aitken <- function(eps = 1e-2) {
  function(model) {
    if (model$ii < 3) return(F)
    #print(model$ii)
    three <- model$logl[c(-2, -1, 0) + model$ii]
    a <- (three[3] - three[2]) / (three[2] - three[1])
    logl_asymp <- three[2] + (three[3] - three[2]) / (1 - a)
    #print(three)
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




EM_control <- function(init_label, criterion, convergence, iter = NULL, iter_max = 500) {
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
    #print(model$ii)
    model$wt  <- Estep(model)
    model     <- Mstep(model)
    model$logl[model$ii] <- Estep(model, logl.only = T)
    #print(convergence(model))
    if (convergence(model)) {
      break
    } else {
      model$ii <- model$ii + 1
    }
    
    
  }
  
  model$logl      <- model$logl[1:min(model$ii, iter_max)]
  model$criterion <- criterion(model)
  #model$converged <- convergence(model)
  #remove(model$ii)
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
  
  model$logl <- logl
  model$criterion <- control$criterion(model)
  model
}


EM <- function(data, G = 1:6, module, control, ...) {
  
  source(module, local = T)
  init_label <- control$init_label
  
  
  all_models <- lapply(G, function(g) {
    
    model <- list()
    model$data  <- data
    model$label <- init_label(data, g)
    model$prop  <- table(model$label) / length(model$label)
    model$wt    <- sapply(1:g, function(ii) model$label == ii) * 1
    model$param <- R.utils::doCall(init_param, ...)(data, model$wt)
    
    if (is.null(control$iter)) {
      model <- EM_until_converge(model, Estep, Mstep, control)
    } else {
      model <- EM_fixed_iter(model, Estep, Mstep, control)
    }
  })
  
  all_models[[which.max(sapply(all_models, function(x) x$criterion))]]
}
