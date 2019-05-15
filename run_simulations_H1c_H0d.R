
source("utility_functions_kdd.R")

nSim <- 1000
nfeat <- 3
parallel <- TRUE
saveCheckPoints <- seq(100, 1000, by = 100)

###################
## weak effect
###################

outFileName <- "sim_study_H1c_H0d_weak.RData"

set.seed(12345)
mySeeds <- sample(seq(1e+4, 1e+5, by = 1), nSim, replace = FALSE)

out <- matrix(NA, nSim, 4)
colnames(out) <- c("AUC.0", 
                   "AUC.u", 
                   "confoundingPval",
                   "cov(C,Y)")

nseq <- seq(300, 500, by = 1)

for (i in seq(nSim)) {
  cat(i, "\n")
  ## simulate data
  set.seed(mySeeds[i])
  p11 <- runif(1, 0.05, 0.45)
  p00 <- p11
  p10 <- 0.5 - p11
  p01 <- 0.5 - p11
  my.n <- sample(nseq, 1)
  my.rho <- runif(1, 0.2, 0.8)
  my.theta <- runif(1, 0.5, 1.0)
  dat <- GenerateData(n = my.n,
                      nfeat, 
                      p11 = p11,
                      p10 = p10,
                      p01 = p01,
                      p00 = p00,
                      alpha = c(0, my.theta), ## beta (label effect), theta (conf effect)
                      rho = my.rho,
                      binVarNames = c("disease", "gender"))
  
  dat$disease <- factor(dat$disease, labels = c("control", "case"))
  dat$gender <- factor(dat$gender, labels = c("female", "male"))
  
  idxTrain <- seq(1, round(my.n/2), by = 1)
  idxTest <- setdiff(seq(my.n), idxTrain)
  ntest <- length(idxTest)
  
  obs <- GetAUC(dat, 
                idxTrain,
                idxTest,
                labelName = "disease",
                featNames = colnames(dat)[-c(1, 2)],
                negClassName = "control", 
                posClassName = "case")
  AUC.0 <- obs$aucObs
  
  aux1 <- RestrictedPermAUC(dat, 
                            idxTrain,
                            idxTest,
                            nperm = ntest,
                            labelName = "disease", 
                            confName = "gender",
                            featNames = colnames(dat)[-c(1, 2)],
                            negClassName = "control", 
                            posClassName = "case",
                            verbose = FALSE,
                            parallel = parallel)
  
  rnull <- aux1$restrictedPermNull
  meanRestrictedPermNull <- mean(rnull)
  varRestrictedPermNull <- var(rnull)
  avarStandNull <- as.numeric(obs$approxVar["v"])
  out[i, "AUC.0"] <- AUC.0
  out[i, "AUC.u"] <- (AUC.0 - meanRestrictedPermNull) * sqrt(avarStandNull/varRestrictedPermNull) + 0.5
  out[i, "confoundingPval"] <- pnorm(meanRestrictedPermNull, 0.5, sqrt(avarStandNull/ntest), lower.tail = FALSE)
  out[i, "cov(C,Y)"] <- p11 * p00 - p10 * p01
  
  if (i %in% saveCheckPoints) {
    cat("saving ", i, " H1 confounding, H0 disease weak simulations", "\n")
    save(out, file = outFileName, compress = TRUE)
  }
}



###################
## moderate effect
###################

outFileName <- "sim_study_H1c_H0d_moderate.RData"

set.seed(12345)
mySeeds <- sample(seq(1e+4, 1e+5, by = 1), nSim, replace = FALSE)

out <- matrix(NA, nSim, 4)
colnames(out) <- c("AUC.0", 
                   "AUC.u", 
                   "confoundingPval",
                   "cov(C,Y)")

nseq <- seq(300, 500, by = 1)

for (i in seq(nSim)) {
  cat(i, "\n")
  ## simulate data
  set.seed(mySeeds[i])
  p11 <- runif(1, 0.05, 0.45)
  p00 <- p11
  p10 <- 0.5 - p11
  p01 <- 0.5 - p11
  my.n <- sample(nseq, 1)
  my.rho <- runif(1, 0.2, 0.8)
  my.theta <- runif(1, 1.0, 1.5)
  dat <- GenerateData(n = my.n,
                      nfeat, 
                      p11 = p11,
                      p10 = p10,
                      p01 = p01,
                      p00 = p00,
                      alpha = c(0, my.theta), ## beta (label effect), theta (conf effect)
                      rho = my.rho,
                      binVarNames = c("disease", "gender"))
  
  dat$disease <- factor(dat$disease, labels = c("control", "case"))
  dat$gender <- factor(dat$gender, labels = c("female", "male"))
  
  idxTrain <- seq(1, round(my.n/2), by = 1)
  idxTest <- setdiff(seq(my.n), idxTrain)
  ntest <- length(idxTest)
  
  obs <- GetAUC(dat, 
                idxTrain,
                idxTest,
                labelName = "disease",
                featNames = colnames(dat)[-c(1, 2)],
                negClassName = "control", 
                posClassName = "case")
  AUC.0 <- obs$aucObs
  
  aux1 <- RestrictedPermAUC(dat, 
                            idxTrain,
                            idxTest,
                            nperm = ntest,
                            labelName = "disease", 
                            confName = "gender",
                            featNames = colnames(dat)[-c(1, 2)],
                            negClassName = "control", 
                            posClassName = "case",
                            verbose = FALSE,
                            parallel = parallel)
  
  rnull <- aux1$restrictedPermNull
  meanRestrictedPermNull <- mean(rnull)
  varRestrictedPermNull <- var(rnull)
  avarStandNull <- as.numeric(obs$approxVar["v"])
  out[i, "AUC.0"] <- AUC.0
  out[i, "AUC.u"] <- (AUC.0 - meanRestrictedPermNull) * sqrt(avarStandNull/varRestrictedPermNull) + 0.5
  out[i, "confoundingPval"] <- pnorm(meanRestrictedPermNull, 0.5, sqrt(avarStandNull/ntest), lower.tail = FALSE)
  out[i, "cov(C,Y)"] <- p11 * p00 - p10 * p01
  
  if (i %in% saveCheckPoints) {
    cat("saving ", i, " H1 confounding, H0 disease moderate simulations", "\n")
    save(out, file = outFileName, compress = TRUE)
  }
}


###################
## moderate effect
###################

outFileName <- "sim_study_H1c_H0d_strong.RData"

set.seed(12345)
mySeeds <- sample(seq(1e+4, 1e+5, by = 1), nSim, replace = FALSE)

out <- matrix(NA, nSim, 4)
colnames(out) <- c("AUC.0", 
                   "AUC.u", 
                   "confoundingPval",
                   "cov(C,Y)")

nseq <- seq(300, 500, by = 1)

for (i in seq(nSim)) {
  cat(i, "\n")
  ## simulate data
  set.seed(mySeeds[i])
  p11 <- runif(1, 0.05, 0.45)
  p00 <- p11
  p10 <- 0.5 - p11
  p01 <- 0.5 - p11
  my.n <- sample(nseq, 1)
  my.rho <- runif(1, 0.2, 0.8)
  my.theta <- runif(1, 1.5, 2)
  dat <- GenerateData(n = my.n,
                      nfeat, 
                      p11 = p11,
                      p10 = p10,
                      p01 = p01,
                      p00 = p00,
                      alpha = c(0, my.theta), ## beta (label effect), theta (conf effect)
                      rho = my.rho,
                      binVarNames = c("disease", "gender"))
  
  dat$disease <- factor(dat$disease, labels = c("control", "case"))
  dat$gender <- factor(dat$gender, labels = c("female", "male"))
  
  idxTrain <- seq(1, round(my.n/2), by = 1)
  idxTest <- setdiff(seq(my.n), idxTrain)
  ntest <- length(idxTest)
  
  obs <- GetAUC(dat, 
                idxTrain,
                idxTest,
                labelName = "disease",
                featNames = colnames(dat)[-c(1, 2)],
                negClassName = "control", 
                posClassName = "case")
  AUC.0 <- obs$aucObs
  
  aux1 <- RestrictedPermAUC(dat, 
                            idxTrain,
                            idxTest,
                            nperm = ntest,
                            labelName = "disease", 
                            confName = "gender",
                            featNames = colnames(dat)[-c(1, 2)],
                            negClassName = "control", 
                            posClassName = "case",
                            verbose = FALSE,
                            parallel = parallel)
  
  rnull <- aux1$restrictedPermNull
  meanRestrictedPermNull <- mean(rnull)
  varRestrictedPermNull <- var(rnull)
  avarStandNull <- as.numeric(obs$approxVar["v"])
  out[i, "AUC.0"] <- AUC.0
  out[i, "AUC.u"] <- (AUC.0 - meanRestrictedPermNull) * sqrt(avarStandNull/varRestrictedPermNull) + 0.5
  out[i, "confoundingPval"] <- pnorm(meanRestrictedPermNull, 0.5, sqrt(avarStandNull/ntest), lower.tail = FALSE)
  out[i, "cov(C,Y)"] <- p11 * p00 - p10 * p01
  
  if (i %in% saveCheckPoints) {
    cat("saving ", i, " H1 confounding, H0 disease strong simulations", "\n")
    save(out, file = outFileName, compress = TRUE)
  }
}

