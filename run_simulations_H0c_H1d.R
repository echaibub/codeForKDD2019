
source("utility_functions_kdd.R")

nSim <- 3000
nfeat <- 3
parallel <- TRUE
saveCheckPoints <- seq(100, 3000, by = 100)

outFileName <- "sim_study_H0c_H1d.RData"

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
  p10 <- p11
  p00 <- 0.5 - p11
  p01 <- 0.5 - p11
  my.n <- sample(nseq, 1)
  my.rho <- runif(1, 0.2, 0.8)
  my.beta <- runif(1, 0.1, 1)
  dat <- GenerateData(n = my.n,
                      nfeat, 
                      p11 = p11,
                      p10 = p10,
                      p01 = p01,
                      p00 = p00,
                      alpha = c(my.beta, 0), ## beta (label effect), theta (conf effect)
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
    cat("saving ", i, " H0 confounding, H1 disease simulations", "\n")
    save(out, file = outFileName, compress = TRUE)
  }
}

