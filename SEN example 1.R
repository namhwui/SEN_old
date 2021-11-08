library(pgmm)
library(qgraph)
source("bareboneEM.R")
data(olive, package = "pgmm")
data <- scale(olive[, -c(1,2)])
label <- olive$Region


set.seed(1)
control <- EM_control(kmeans_label, criterion_BIC, convergence_aitken(eps = 1e-2), iter_max = 1000)
model1 <- EM(data, 3, "SEN functions.R", control, alpha = 0.5, rho = 0, q = 1)
model2 <- EM(data, 3, "SEN functions.R", control, alpha = 0.9, rho = 1, q = 1)
model3 <- EM(data, 3, "SEN functions.R", control, alpha = 0.5, rho = 1, q = 1)
model4 <- EM(data, 3, "SEN functions.R", control, alpha = 0, rho = 1, q = 1)


# log-likelihood values
c(tail(model1$logl, 1), tail(model2$logl, 1), tail(model3$logl, 1), tail(model4$logl, 1))

# Adjusted Rand Index
# install e1071 package first
c(e1071::classAgreement(table(label, model1$label))$crand,
  e1071::classAgreement(table(label, model2$label))$crand,
  e1071::classAgreement(table(label, model3$label))$crand,
  e1071::classAgreement(table(label, model4$label))$crand)


# correlation graph
# install qgraph package first
par(mfcol = c(3,4))
qgraph(cov2cor(model1$param[[1]]$sigma), minimum = 0.5)
qgraph(cov2cor(model1$param[[2]]$sigma), minimum = 0.5)
qgraph(cov2cor(model1$param[[3]]$sigma), minimum = 0.5)

qgraph(cov2cor(model2$param[[1]]$sigma), minimum = 0.5)
qgraph(cov2cor(model2$param[[2]]$sigma), minimum = 0.5)
qgraph(cov2cor(model2$param[[3]]$sigma), minimum = 0.5)

qgraph(cov2cor(model3$param[[1]]$sigma), minimum = 0.5)
qgraph(cov2cor(model3$param[[2]]$sigma), minimum = 0.5)
qgraph(cov2cor(model3$param[[3]]$sigma), minimum = 0.5)

qgraph(cov2cor(model4$param[[1]]$sigma), minimum = 0.5)
qgraph(cov2cor(model4$param[[2]]$sigma), minimum = 0.5)
qgraph(cov2cor(model4$param[[3]]$sigma), minimum = 0.5)




