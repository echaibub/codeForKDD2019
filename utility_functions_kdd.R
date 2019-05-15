
## Load dependencies.
##
library(MASS)
library(dplyr)
library(plyr)
library(pROC)
library(randomForest)
library(doMC)
registerDoMC(detectCores() - 2)


######################################
## Functions used to simulate 
## the synthetic data
###################################### 

## Creates the correlation matrix
## used to generate correlated 
## features.
##
CreateCorrelationMatrix <- function(rho, p) {
  aux1 <- matrix(rep(1:p, p), p, p)
  aux2 <- matrix(rep(1:p, each = p), p, p) 
  
  rho^abs(aux1 - aux2)
}

## Draw samples from the bivariate
## Bernoulli distribution.
##
SampleBivariateBernoulli <- function(n, p00, p01, p10, p11, binVarNames = NULL) {
  aux <- sample(c("0,0", "0,1", "1,0", "1,1"), n, replace = TRUE, prob = c(p00, p01, p10, p11))
  out <- matrix(NA, n, 2)
  idx00 <- which(aux == "0,0")
  out[idx00, 1] <- 0
  out[idx00, 2] <- 0
  idx01 <- which(aux == "0,1")
  out[idx01, 1] <- 0
  out[idx01, 2] <- 1
  idx10 <- which(aux == "1,0")
  out[idx10, 1] <- 1
  out[idx10, 2] <- 0
  idx11 <- which(aux == "1,1")
  out[idx11, 1] <- 1
  out[idx11, 2] <- 1
  colnames(out) <- binVarNames
  
  out
}

## Simulate the synthetic data.
##
GenerateData <- function(n,
                         nfeat, 
                         p11,
                         p10,
                         p01,
                         p00,
                         alpha, ## c(beta, theta) [disease, confounding] effects
                         rho = 0.5,
                         binVarNames = NULL) {
  X <- SampleBivariateBernoulli(n, p00, p01, p10, p11, binVarNames)
  mu <- rep(0, nfeat)
  Sig <- CreateCorrelationMatrix(rho, p = nfeat)
  err <- mvrnorm(n, mu, Sig)
  lin_pred <- X %*% alpha
  f <- lin_pred[, 1] + err
  colnames(f) <- paste("f", seq(nfeat), sep = "")
  
  data.frame(X, f)
}


#########################################
## Functions used in the analyses
#########################################

## Computes observed AUC (from a random forest fit).
## Also returns the variance of the 
## normal approximation for the standard permutation
## null distribution.
##
GetAUC <- function(dat,
                   idxTrain, 
                   idxTest, 
                   labelName, 
                   featNames,
                   negClassName, 
                   posClassName) {
  dat <- dat[, c(labelName, featNames)]
  dat[, labelName] <- factor(as.character(dat[, labelName]), 
                             levels = c(negClassName, posClassName)) 
  myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
  fit <- randomForest(myFormula, data = dat[idxTrain,], ntree = 1000)
  predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "prob")
  rocObj <- roc(dat[idxTest, 1], predProbs[, posClassName], direction = "<", 
                levels = c(negClassName, posClassName))    
  aucObs <- pROC::auc(rocObj)[1]
  approxVar <- GetNormApproxVarAUC(dat[idxTest, 1], negClassName, posClassName)
  
  list(aucObs = aucObs, predProbs = predProbs, rocObj = rocObj, approxVar = approxVar)
}

## Calculates the variance of the normal approximation 
## for the standard permutation null distribution.
##
GetNormApproxVarAUC <- function(ytest, negClassName, posClassName) {
  ytest <- factor(ytest)
  ylevels <- levels(ytest)
  n1 <- sum(ytest == negClassName)
  n2 <- sum(ytest == posClassName)
  n <- n1 + n2
  v <- (n + 1)/(12 * n1 * n2)
  
  c(v = v, n = n, nNeg = n1, nPos = n2)
}

## Performs one restricted permutation of
## the data.
##
RestrictedResponseShuffling <- function(dat, respName, confName) {
  sdat <- dat
  confLevels <- unique(dat[, confName])
  nlevels <- length(confLevels)
  for (i in seq(nlevels)) {
    idx <- which(dat[, confName] == confLevels[i])
    shuffledLabels <- dat[idx, respName][sample(length(idx))]
    sdat[idx, respName] <- shuffledLabels
  }
  
  sdat
}

## Generates the restricted permutation
## null distribution for the AUC metric
## for the random forest classifier.
##
RestrictedPermAUC <- function(dat, 
                             idxTrain,
                             idxTest,
                             nperm,
                             labelName, 
                             confName,
                             featNames,
                             negClassName, 
                             posClassName,
                             verbose = FALSE,
                             parallel = TRUE) {
  dat <- dat[, c(labelName, confName, featNames)]
  dat[, labelName] <- factor(as.character(dat[, labelName]), 
                             levels = c(negClassName, posClassName)) 
  trainData <- dat[idxTrain,]
  testData <- dat[idxTest,]
  myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
  if (verbose) {
    res_auc <- plyr::llply(1:nperm, .parallel = parallel,  function(num){
      trainDataS <- RestrictedResponseShuffling(trainData, labelName, confName)
      testDataS <- RestrictedResponseShuffling(testData, labelName, confName)
      fitS <- randomForest(myFormula, data = trainDataS)
      predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "prob")
      rocObjS <- roc(testDataS[, 1], predProbsS[, posClassName], direction = "<", 
                     levels = c(negClassName, posClassName)) 
      pROC::auc(rocObjS)[1]
    }, .progress = progress_text(char = "."))
  }
  else {
    res_auc <- plyr::llply(1:nperm, .parallel = parallel,  function(num){
      trainDataS <- RestrictedResponseShuffling(trainData, labelName, confName)
      testDataS <- RestrictedResponseShuffling(testData, labelName, confName)
      fitS <- randomForest(myFormula, data = trainDataS)
      predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "prob")
      rocObjS <- roc(testDataS[, 1], predProbsS[, posClassName], direction = "<", 
                     levels = c(negClassName, posClassName)) 
      pROC::auc(rocObjS)[1]
    })
  }
  restrictedPermNull <- unlist(res_auc)
  
  list(restrictedPermNull = restrictedPermNull)
}

## Performs a random split of the data into training and test sets.
##
TrainTestStratifiedSplit <- function(dat, ntest, respName, confName) {
  tb <- as.data.frame.matrix(table(dat[, confName], dat[, respName]))
  n <- nrow(dat)
  ntrain <- n - ntest
  probs <- tb/n
  testNumbers <- round(probs * ntest)
  trainNumbers <- round(probs * ntrain)
  confLevels <- rownames(tb)
  respLevels <- colnames(tb)
  idxTest <- c()
  for (i in seq(length(confLevels))) {
    for (j in seq(length(respLevels))) {
      aux <- which(dat[, confName] == confLevels[i] & dat[, respName] == respLevels[j])
      idxTest <- c(idxTest, sample(aux, testNumbers[i, j], replace = FALSE))
    }
  }
  idxTrain <- setdiff(seq(n), idxTest)
  
  list(idxTrain = idxTrain, 
       idxTest = idxTest, 
       probs = probs, 
       trainNumbers = trainNumbers,
       testNumbers = testNumbers)
}

## Generates the baseline train and test sets
## for the synthetic data illustrations presented 
## in the "Accounting for the confounder/response 
## association structure in the target population"
## section.
##
## freqs is a matrix
##
GetBaselineTrainTestSets <- function(dat, freqs, respName, confName) {
  GetTestIndexes <- function(dat, freqs, respName, confName) {
    confLevels <- rownames(freqs)
    respLevels <- colnames(freqs)
    idxTest <- c()
    for (i in seq(length(confLevels))) {
      for (j in seq(length(respLevels))) {
        aux <- which(dat[, confName] == confLevels[i] & dat[, respName] == respLevels[j])
        idxTest <- c(idxTest, sample(aux, freqs[i, j], replace = FALSE))
      }
    }
    idxTest
  }
  idxTest <- GetTestIndexes(dat, freqs, respName = "disease", confName = "gender")
  datTest <- dat[idxTest,]
  datTrain <- dat[-idxTest,]
  idxTrain <- GetTestIndexes(datTrain, freqs, respName = "disease", confName = "gender")
  datTrain <- datTrain[idxTrain,]
  
  list(datTrain = datTrain, datTest = datTest)
}








