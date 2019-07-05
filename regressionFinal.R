# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to work with the output plink and admixture on biolync.case.edu server. The goal is to compare two implementations of lasso regression. One with 

require(tidyverse)
require(ssh)
require(glmnet)
require(KRIS)
require(glmnet)
library(GGally)

# Connect to biolync server
session <- ssh_connect("jxl2059@biolync.case.edu")

## Load Admixture results
load(file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/admixture2")
load(file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/ethFam")
## Load initial lasso results
load(file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/plinkPedJoin")

# Join with admixture results
obs <- merge(plinkPedJoin, select(admixture2,  IID, CEU, YRI))
# Perform multivariate linear regression with ridge regression vs lasso
# Generate observation, obs, and response, resp, matrices
# Filter for only relevant snps rs1052406 and rs1009848
obs <- select(obs, phenotype, rs1052406, rs1009848, AGE, SEX, CEU, YRI)
pairData <- obs
# obs <- data.matrix(obs)
obs <- model.matrix(phenotype ~ ., obs)[,-1]
resp <- data.matrix(select(plinkPedJoin, phenotype))

# lambda values to test
grid <- 10^seq(10, -2, length = 100)

# Lasso-Admixture 
# Lasso, alpha = 1, with Admixture
lasso.fit <- glmnet(obs, resp, alpha = 1, lambda = grid)
summary(lasso.fit)
# see coefficients with given lambda
# coef(lasso.fit)[,100]
jpeg('figures/lassoAdmixtureCoeff.jpeg')
plot(lasso.fit, label = TRUE)
dev.off()
set.seed(1)
# LOOCV since nfolds = number of observations
cv.lasso <- cv.glmnet(obs, resp, alpha = 1, nfolds = 15, grouped = FALSE)
jpeg('figures/lassoAdmixtureCV.jpeg')
plot(cv.lasso)
dev.off()
bestLambda <-  cv.lasso$lambda.min
bestLambda
# 0.01396051
predict(lasso.fit, type = "coefficients", s = bestLambda)
# rs1052406       -0.04389696
lasso.fit.best <- glmnet(obs, resp, alpha = 1, lambda = bestLambda)
# Evaluation
# MSE
# [1] 0.002557937
lasso.pred <- predict(lasso.fit ,s = bestLambda, newx = obs)
mean((lasso.pred - resp)^2)
# AIC
tLL <- lasso.fit.best$nulldev - deviance(lasso.fit.best)
k <- lasso.fit.best$df
n <- lasso.fit.best$nobs
AICc <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
AICc
# [1] 2.283352

### PCA
# Load PCE eignenvectors
load(file = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/eigenMerge")
eigenMerge$IID = row.names(eigenMerge)
eigenMerge <- select(eigenMerge, IID, PC1, PC2)
# Join with PCA results
obs <- merge(plinkPedJoin, select(eigenMerge,  IID, PC1, PC2), by = "IID")
resp <- obs$phenotype
# Generate observation, obs, and response, resp, matrices
# Filter for only relevant snps rs1052406 and rs1009848
obs <- select(obs, phenotype, rs1052406, rs1009848, AGE, SEX, PC1, PC2)
# obs <- data.matrix(obs)
obs <- model.matrix(phenotype ~ ., obs)[,-1]
resp <- data.matrix(select(plinkPedJoin, phenotype))


# Lasso-PCA 
# Lasso, alpha = 1, with PCA
lasso.fit <- glmnet(obs, resp, alpha = 1, lambda = grid)
summary(lasso.fit)
# see coefficients with given lambda
# coef(lasso.fit)[,100]
jpeg('figures/lassoPcaCoeff.jpeg')
plot(lasso.fit, label = TRUE)
dev.off()
set.seed(1)
# LOOCV since nfolds = number of observations
cv.lasso <- cv.glmnet(obs, resp, alpha = 1, nfolds = 15, grouped = FALSE)
jpeg('figures/lassoPcaCV.jpeg')
plot(cv.lasso)
dev.off()
bestLambda <-  cv.lasso$lambda.min
bestLambda
#  0.013960511
predict(lasso.fit, type = "coefficients", s = bestLambda)
# (Intercept)  0.15516923
# rs1052406   -0.04379660
# rs1009848    .         
# AGE          .         
# SEX1         .         
# PC1          .         
# PC2         -0.06568802
lasso.fit.best <- glmnet(obs, resp, alpha = 1, lambda = bestLambda)
# Evaluation
# MSE
lasso.pred <- predict(lasso.fit ,s = bestLambda, newx = obs)
mean((lasso.pred - resp)^2)
# [1] 0.00251772
# AIC
tLL <- lasso.fit.best$nulldev - deviance(lasso.fit.best)
k <- lasso.fit.best$df
n <- lasso.fit.best$nobs
AICc <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
AICc
# 4.975214

# Create pairs plot
pairData <- cbind(pairData, as.data.frame(obs)$PC1, as.data.frame(obs)$PC2)
colnames(pairData)[8:9] <- c("PC1", "PC2")
pairData$rs1052406 <- as.factor(as.character(pairData$rs1052406))
jpeg('figures/pairsPlot.jpeg')
ggpairs(pairData, ggplot2::aes(colour = rs1052406))
dev.off()

