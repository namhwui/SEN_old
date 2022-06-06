# Clustering the Olive data 
data(olive)
dat <- scale(olive[, -c(1, 2)])
region <- olive$Region
area <- olive$Area

set.seed(1)
# run until convergence
# k-means initialization
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 1, rho = 10)
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)

# run until convergence
# hierarchical initialization
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 1, rho = 10, 
                control = EM_control(init_label = init_hc))
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)

# fixed number of iterations
# k-means initialization
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 1, rho = 10,
                control = EM_control(iter = 200))
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)

# fixed number of iterations
# k-means initialization
# L2 penalty only
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 0, rho = 10)
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)


# fixed number of iterations
# k-means initialization
# Both L1 and L2 penalties
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 0.5, rho = 10)
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)


# fixed number of iterations
# k-means initialization
# weaker penalty
model <- SEN_EM(data, G = 1:4, q = 4, alpha = 1, rho = 2)
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)


# fixed number of iterations
# k-means initialization
# fewer factors
model <- SEN_EM(data, G = 1:4, q = 2, alpha = 1, rho = 10)
ARI(region, model$label, print_table = T)
ARI(area, model$label, print_table = T)
