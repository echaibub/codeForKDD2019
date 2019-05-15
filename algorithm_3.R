
source("utility_functions_kdd.R")

##################################
## Additional utility functions.
##################################

ConfounderShuffling <- function(dat, confName) {
  n <- nrow(dat)
  dat[, confName] <- dat[sample(n), confName]
  
  dat
}

ConfPermTestAUC <- function(dat, 
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
  confPermNull <- rep(NA, nperm)
  ## get observed AUC
  obs <- GetAUC(dat, 
                idxTrain,
                idxTest,
                labelName,
                featNames,
                negClassName, 
                posClassName)
  approxVar <- obs$approxVar
  AUC.0 <- obs$aucObs
  
  trainData <- dat[idxTrain,]
  testData <- dat[idxTest,]
  
  npermr <- nrow(testData)
  myFormula <- as.formula(paste(labelName, " ~ ", paste(featNames, collapse = " + ")))
  for (i in seq(nperm)) {
    cat(i, "\n")
    trainDataP <- ConfounderShuffling(trainData, confName)
    testDataP <- ConfounderShuffling(testData, confName)
    if (verbose) {
      res_auc <- plyr::llply(1:npermr, .parallel = parallel,  function(num){
        trainDataS <- RestrictedResponseShuffling(trainDataP, labelName, confName)
        testDataS <- RestrictedResponseShuffling(testDataP, labelName, confName)
        fitS <- randomForest(myFormula, data = trainDataS)
        predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "prob")
        rocObjS <- roc(testDataS[, 1], predProbsS[, posClassName], direction = "<", 
                       levels = c(negClassName, posClassName)) 
        pROC::auc(rocObjS)[1]
      }, .progress = progress_text(char = "."))
    }
    else {
      res_auc <- plyr::llply(1:npermr, .parallel = parallel,  function(num){
        trainDataS <- RestrictedResponseShuffling(trainDataP, labelName, confName)
        testDataS <- RestrictedResponseShuffling(testDataP, labelName, confName)
        fitS <- randomForest(myFormula, data = trainDataS)
        predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "prob")
        rocObjS <- roc(testDataS[, 1], predProbsS[, posClassName], direction = "<", 
                       levels = c(negClassName, posClassName)) 
        pROC::auc(rocObjS)[1]
      })
    }
    restrictedPermNull <- unlist(res_auc)
    confPermNull[i] <- mean(restrictedPermNull)
  }
  
  list(AUC.0 = AUC.0,
       npermr = npermr,
       approxVar = approxVar,
       confPermNull = confPermNull)
}


###############################################
###############################################
###############################################

myseed <- 123456789

set.seed(myseed)
dat1 <- GenerateData(n = 20,
                     nfeat = 5, 
                     p11 = 0.40,
                     p10 = 0.15,
                     p01 = 0.15,
                     p00 = 0.30,
                     alpha = c(1, 1),
                     rho = 0.5,
                     binVarNames = c("disease", "gender"))
dat1$disease <- factor(dat1$disease, labels = c("control", "case"))
dat1$gender <- factor(dat1$gender, labels = c("female", "male"))

table(dat1[1:10,]$gender, dat1[1:10,]$disease)
table(dat1[11:20,]$gender, dat1[11:20,]$disease)



set.seed(myseed)
dat2 <- GenerateData(n = 200,
                     nfeat = 5, 
                     p11 = 0.40,
                     p10 = 0.15,
                     p01 = 0.15,
                     p00 = 0.30,
                     alpha = c(1, 1),
                     rho = 0.5,
                     binVarNames = c("disease", "gender"))
dat2$disease <- factor(dat2$disease, labels = c("control", "case"))
dat2$gender <- factor(dat2$gender, labels = c("female", "male"))


set.seed(myseed)
aux1 <- ConfPermTestAUC(dat1, 
                        idxTrain = seq(1, 10, by = 1),
                        idxTest = seq(11, 20, by = 1),
                        nperm = 1e+4,
                        labelName = "disease", 
                        confName = "gender",
                        featNames = paste("f", 1:5, sep = ""),
                        negClassName = "control", 
                        posClassName = "case",
                        verbose = FALSE,
                        parallel = TRUE)


set.seed(myseed)
aux2 <- ConfPermTestAUC(dat2, 
                        idxTrain = seq(1, 100, by = 1),
                        idxTest = seq(101, 200, by = 1),
                        nperm = 1e+4,
                        labelName = "disease", 
                        confName = "gender",
                        featNames = paste("f", 1:5, sep = ""),
                        negClassName = "control", 
                        posClassName = "case",
                        verbose = FALSE,
                        parallel = TRUE)

#save(aux1, aux2, file = "outputs_algo3.RData", compress = TRUE)

xaxis1 <- seq(min(aux1$confPermNull), max(aux1$confPermNull), length.out = 1000)
adensi1 <- dnorm(xaxis1, 0.5, sqrt(aux1$approxVar["v"]/aux1$npermr))

xaxis2 <- seq(min(aux2$confPermNull), max(aux2$confPermNull), length.out = 1000)
adensi2 <- dnorm(xaxis2, 0.5, sqrt(aux2$approxVar["v"]/aux2$npermr))

densi1 <- density(aux1$confPermNull)
densi2 <- density(aux2$confPermNull)


##################################
## Plot Supplementary Figure S1. 
##################################

load("outputs_algo3.RData")

xaxis1 <- seq(min(aux1$confPermNull), max(aux1$confPermNull), length.out = 1000)
adensi1 <- dnorm(xaxis1, 0.5, sqrt(aux1$approxVar["v"]/aux1$npermr))

xaxis2 <- seq(min(aux2$confPermNull), max(aux2$confPermNull), length.out = 1000)
adensi2 <- dnorm(xaxis2, 0.5, sqrt(aux2$approxVar["v"]/aux2$npermr))

densi1 <- density(aux1$confPermNull)
densi2 <- density(aux2$confPermNull)

par(mfrow = c(2, 1), mar = c(3, 2.5, 1.5, 0.5), mgp = c(1.5, 0.5, 0))
plot(densi1, xlim = c(0.15, 0.85), ylim = c(0, 6.1), col = "blue",
     xlab = "average AUC", main = "Test set size = 10")
lines(xaxis1, adensi1, col = "red")
legend("topleft", bty = "n", text.col = c("blue", "red"),
       legend = c("permutation null distribution", "normal approximation null"))
####
plot(densi2, xlim = c(0.47, 0.53), col = "blue",
     xlab = "average AUC", main = "Test set size = 100")
lines(xaxis2, adensi2, col = "red")
legend("topleft", bty = "n", text.col = c("blue", "red"),
       legend = c("permutation null distribution", "normal approximation null"))
par(mfrow = c(1, 1))

