# Copyright: CC-BY-SA 4.0
# Author Notes:  Please note this work is in conjunction with Dr. Dana Crawford and Dr. William Bush at Case Western Reserve Universithy and Cleveland Institute of Computational Biology. The data used below is not my own.
# File Description: The code below is intended to work with the output plink on biolync.case.edu server. The goal is to create a multivariate model using age, sex, and genotypes as covariates to predict the response, productive clonality of T Cells.

require(tidyverse)
require(ssh)
require(glmnet)
require(KRIS)

# Connect to biolync server
session <- ssh_connect("jxl2059@biolync.case.edu")

# Read in  plink5 results
plinkBed <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink5/window.bed", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkBim <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink5/window.bim", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkFam <- scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/plinkResults/plinkFiltering/plink5/window.fam", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
plinkRead <- read.bed("/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/window.bed", "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/window.bim", "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/window.fam", 
                      only.snp = FALSE)
# Read in tcrEmrPheno (phenotype information)
scp_download(session, "/storage/mips/MIPS_Updated.2019-02-21/jxl2059/phenotypes/tcrEmrPheno.txt", to = "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data")
tcrEmrPheno <- read.delim(
  "/Users/linjo/GoogleDrive/CaseWesternUniversity/tcr-project/data/tcrEmrPheno.txt", 
  sep = "\t", quote = "")

# Creating separate data frames from bed object
plinkPed <- plinkRead$snp
plinkPed <- as.data.frame(plinkPed)
plinkMap <- plinkRead$snp.info
plinkFam <- plinkRead$ind.info

# Add SNP column names to plinkPed
# TODO: Verify order of snps from plinkMap can be transposed
colnames(plinkPed) <- plinkMap$ID
# Substitute commas for pipes in column names of snps
colnames(plinkPed)[grep(",", colnames(plinkPed))] <- 
  gsub(",","|",colnames(plinkPed)[grep(",", colnames(plinkPed))])
# Add FID and IID columns
# TODO: Verify we can just bind the plinkFam columns - the order is the same?
plinkPed$FID <- plinkFam$FamID
plinkPed$IID <- plinkFam$IndID
# Also, add in phenotypes to plinkPed
# TODO: Verify - the order is the same?
plinkPed$phenotype <- plinkFam$phenotype


# Create other covariates for AGE and SEX
tcrEmrPhenoAgeSex <- select(tcrEmrPheno, FID, IID, AGE, SEX)

# Merge other covariates with plinkPed
plinkPedJoin <- full_join(plinkPed, tcrEmrPhenoAgeSex, by = c("FID", "IID"))

# Reorder
plinkPedJoin <- select(plinkPedJoin, FID, IID, AGE, SEX, phenotype, everything())

# Recode sex variable as categorical
# Males are 1
plinkPedJoin$SEX[which(plinkPedJoin$SEX == 1)] = 1
# Females are 0
plinkPedJoin$SEX[which(plinkPedJoin$SEX == 2)] = 0
plinkPedJoin$SEX <- as.factor(plinkPedJoin$SEX)

# Perform multivariate linear regression with ridge regression vs lasso
# Generate observation, obs, and response, resp, matrices
obs <- select(plinkPedJoin, -FID, -IID, -phenotype, -AGE, -SEX)
# Filter out columns where there are NA values
which(colSums(is.na(obs)) != 0)
# these ones are removed JHU_7.142361724, rs151232837, rs6947539 
obs <- obs[,colSums(is.na(obs)) == 0]
obs <- data.matrix(obs)
resp <- data.matrix(select(plinkPedJoin, phenotype))

# lambda values to test
grid <- 10^seq(10, -2, length = 100)

# Ridge regression, alpha = 0
ridge.fit <- glmnet(obs, resp, alpha = 0, lambda = grid)
summary(ridge.fit)
# see coefficients with given lambda
# coef(ridge.fit)[,100]
jpeg('figures/ridgeCoeff1.jpeg')
plot(ridge.fit, label = TRUE)
dev.off()
set.seed(1)
# LOOCV since nfolds = number of observations
cv.ridge <- cv.glmnet(obs, resp, alpha = 0, nfolds = 15, grouped = FALSE)
jpeg('figures/ridge1.jpeg')
plot(cv.ridge)
dev.off()
bestLambda <-  cv.ridge$lambda.min
bestLambda
# 29.68392
predict(ridge.fit, type = "coefficients", s = bestLambda)

# Lasso, alpha = 1
lasso.fit <- glmnet(obs, resp, alpha = 1, lambda = grid)
summary(lasso.fit)
# see coefficients with given lambda
# coef(lasso.fit)[,100]
jpeg('figures/lassoCoeff1.jpeg')
plot(lasso.fit, label = TRUE)
dev.off()
set.seed(1)
# LOOCV since nfolds = number of observations
cv.lasso <- cv.glmnet(obs, resp, alpha = 1, nfolds = 15, grouped = FALSE)
jpeg('figures/lasso1.jpeg')
plot(cv.lasso)
dev.off()
bestLambda <-  cv.lasso$lambda.min
bestLambda
# 0.03575439
predict(lasso.fit, type = "coefficients", s = bestLambda)
# rs1052406       -0.002969940
# rs1009848        0.002077136



