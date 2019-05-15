
source("utility_functions_kdd.R")

####################################
## Additional utility functions
####################################

## Convert the longitudinal feature measurements
## of each participant to a single value 
## (the median across all longitudinal values).
##
CollapseFeaturesMedian <- function(x, labelName, covNames, subjectIdName, featNames) {
  ids <- as.character(unique(x[, subjectIdName]))
  nids <- length(ids)
  nfeat <- length(featNames)
  nvar <- length(c(subjectIdName, labelName, covNames))
  out <- data.frame(matrix(NA, nids, nfeat + nvar))
  colnames(out) <- c(subjectIdName, labelName, covNames, featNames)
  rownames(out) <- ids
  for (i in seq(nids)) {
    sdat <- x[which(x[, subjectIdName] == ids[i]),]
    out[i, subjectIdName] <- ids[i]
    out[i, labelName] <- as.character(sdat[1, labelName])
    out[i, covNames] <- sdat[1, covNames]
    out[i, (nvar + 1):(nfeat + nvar)] <- apply(sdat[, featNames], 2, median, na.rm = TRUE)
  }
  
  out
}

## Get the demographic data for each participant.
## (take a single value from a matrix containing
## duplicated values for each participant).
##
CollapseDemo <- function(x, subjectIdName) {
  ids <- as.character(unique(x[, subjectIdName]))
  nids <- length(ids)
  nc <- ncol(x)
  out <- data.frame(matrix(NA, nids, nc))
  colnames(out) <- names(x)
  rownames(out) <- ids
  for (i in seq(nids)) {
    sdat <- x[which(x[, subjectIdName] == ids[i]),]
    out[i,] <- sdat[1,]
  }
  
  out
} 

## Discretize the age variable into nlevels discrete levels.
##
GetDiscretizedAge <- function(dat, nlevels, levelNames = NULL, breaks = NULL) {
  if (is.null(breaks)) {
    if (!is.null(levelNames)) {
      out <- cut(dat$age, breaks = nlevels, labels = levelNames)
    }
    else {
      levelNames <- paste("age", seq(nlevels), sep = "")
      out <- cut(dat$age, breaks = nlevels, labels = levelNames)
    }
  }
  else {
    if (!is.null(levelNames)) {
      out <- cut(dat$age, breaks = breaks, labels = levelNames)
    }
    else {
      nlevels <- length(breaks) + 1
      levelNames <- paste("age", seq(nlevels), sep = "")
      out <- cut(dat$age, breaks = breaks, labels = levelNames)
    }
  }
  
  out
}

## Computes observed AUC (from a logistic regression fit).
## Also returns the variance of the 
## normal approximation for the standard permutation
## null distribution.
##
GetAUCGlm <- function(dat,
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
  fit <- glm(myFormula, data = dat[idxTrain,], family = "binomial")
  predProbs <- predict(fit, dat[idxTest, -1, drop = FALSE], type = "response")
  rocObj <- roc(dat[idxTest, 1], predProbs, direction = "<", 
                levels = c(negClassName, posClassName))    
  aucObs <- pROC::auc(rocObj)[1]
  approxVar <- GetNormApproxVarAUC(dat[idxTest, 1], negClassName, posClassName)
  
  list(aucObs = aucObs, predProbs = predProbs, rocObj = rocObj, approxVar = approxVar)
}

## Generates the restricted permutation
## null distribution for the AUC metric
## for the logistic regression classifier.
##
RestrictedPermAUCGlm <- function(dat, 
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
      fitS <- glm(myFormula, data = trainDataS, family = "binomial")
      predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "response")
      rocObjS <- roc(testDataS[, 1], predProbsS, direction = "<", 
                     levels = c(negClassName, posClassName)) 
      pROC::auc(rocObjS)[1]
    }, .progress = progress_text(char = "."))
  }
  else {
    res_auc <- plyr::llply(1:nperm, .parallel = parallel,  function(num){
      trainDataS <- RestrictedResponseShuffling(trainData, labelName, confName)
      testDataS <- RestrictedResponseShuffling(testData, labelName, confName)
      fitS <- glm(myFormula, data = trainDataS, family = "binomial")
      predProbsS <- predict(fitS, testDataS[, -1, drop = FALSE], type = "response")
      rocObjS <- roc(testDataS[, 1], predProbsS, direction = "<", 
                     levels = c(negClassName, posClassName)) 
      pROC::auc(rocObjS)[1]
    })
  }
  restrictedPermNull <- unlist(res_auc)
  
  list(restrictedPermNull = restrictedPermNull)
}

## Generate the matched data.
##
MatchedDataByConfounder <- function(dat, labelName, confName) {
  labelLevels <- levels(factor(dat[, labelName]))
  confLevels <- levels(factor(dat[, confName]))
  nConfLevels <- length(confLevels)
  idxList1 <- vector(mode = "list", length = nConfLevels)
  idxList2 <- idxList1
  auxLengths1 <- rep(NA, nConfLevels)
  auxLengths2 <- auxLengths1
  for (j in seq(nConfLevels)) {
    idxList1[[j]] <- which(dat[, labelName] == labelLevels[1] & dat[, confName] == confLevels[j])
    auxLengths1[j] <- length(idxList1[[j]])
    idxList2[[j]] <- which(dat[, labelName] == labelLevels[2] & dat[, confName] == confLevels[j])
    auxLengths2[j] <- length(idxList2[[j]])
  }
  n2keep1 <- min(auxLengths1)
  n2keep2 <- min(auxLengths2)
  idx2keep1 <- vector(mode = "list", length = nConfLevels)
  idx2keep2 <- idx2keep1
  for (j in seq(nConfLevels)) {
    idx2keep1[[j]] <- sample(idxList1[[j]], n2keep1, replace = FALSE)
    idx2keep2[[j]] <- sample(idxList2[[j]], n2keep2, replace = FALSE)
  }  
  keep <- sort(c(unlist(idx2keep1), unlist(idx2keep2)))
  
  dat[keep,]
}

## Generates the augmented dataset using the 
## approximate IPW method based in prosensity
## scores computed using logistic regression.
##
IPWPseudoPopulation <- function(dat, covNames) {
  dat$PD <- as.numeric(as.factor(dat$professional.diagnosis)) - 1
  myformula <- as.formula(paste("PD", " ~ ", paste(covNames, collapse = " + ")))
  fit <- glm(myformula, data = dat, family = "binomial")
  probs <- predict(fit, type = "response")
  wi <- 1/probs
  idx0 <- dat$PD == 0
  wi[idx0] <- 1 - probs[idx0]
  ti <- round(1/wi)
  ti[ti == 0] <- 1
  idx <- rep(seq(nrow(dat)), times = ti)
  rdat <- dat[idx,]
  rdat <- rdat[, -which(names(rdat) == "PD")]
  
  list(rdat = rdat, ti = ti, wi = wi)
}

## Computes the unconfounded AUC score.
##
GetUnconfoundedAUC <- function(obs, permNull) {
  auc.0 <- obs$aucObs
  avarStandNull <- as.numeric(obs$approxVar["v"])
  auc.u <- (auc.0 - mean(permNull)) * sqrt(avarStandNull)/sd(permNull) + 0.5
  
  auc.u
}

## Computes the confounding test p-value.
##
ApproximateConfoundingPval <- function(obs, permNull) {
  obsStat <- mean(permNull)
  n <- obs$approxVar["n"]
  nNeg <- obs$approxVar["nNeg"]
  nPos <- obs$approxVar["nPos"]
  v <- (nNeg + nPos + 1)/(12 * nNeg * nPos * n)
  
  pnorm(obsStat, 0.5, sqrt(v), lower.tail = FALSE)
}



#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################

## Load the Synapse client (that hosts
## the mPower data).
##
require(synapseClient)
synapseLogin()

## Note that in order to get access to the 
## mPower training data a researcher needs to
## first apply for data access. 
##
## Visit
## https://www.synapse.org/#!Synapse:syn4993293/wiki/247859
## for details on how to get access.

load(getFileLocation(synGet("syn12182549")))
feats <- read.csv(getFileLocation(synGet("syn11027900")), header = TRUE)

#####################
## Shape the data.
#####################

dat <- data.frame(mDat, feats[match(mDat$recordId, feats$recordId),])
dat$age <- as.numeric(dat$age)
dat <- na.omit(dat)
idx_train <- which(dat$datSplit == "train")
idx_test <- which(dat$datSplit == "test")
dat_train <- dat[idx_train,]
dat_test <- dat[idx_test,]
featNames <- colnames(feats)[-1]

## Collapse the longitudinal features.
##
train_dat_c <- CollapseFeaturesMedian(x = dat_train, 
                                      labelName = "professional.diagnosis", 
                                      covNames = c("age", "gender"), 
                                      subjectIdName = "healthCode", 
                                      featNames = featNames)
test_dat_c <- CollapseFeaturesMedian(x = dat_test, 
                                     labelName = "professional.diagnosis", 
                                     covNames = c("age", "gender"), 
                                     subjectIdName = "healthCode", 
                                     featNames = featNames)

## Generate the discretized age variable.
##
train_dat_c$dage <- GetDiscretizedAge(train_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                      levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))
test_dat_c$dage <- GetDiscretizedAge(test_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                     levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))

## Create the combined confounder.
##
train_dat_c$combinedConfounder <- paste(train_dat_c$dage, train_dat_c$gender, sep = "_")
test_dat_c$combinedConfounder <- paste(test_dat_c$dage, test_dat_c$gender, sep = "_")

dat_c <- rbind(train_dat_c, test_dat_c)

nperm <- 1e+4

#########################################
## Run the analyses with no adjustment.
#########################################

itrain <- seq(nrow(train_dat_c))
itest <- max(itrain) + seq(nrow(test_dat_c))

## Get the observed AUC score
## using the random forest classifier.
##
set.seed(1234567)
obsRf <- GetAUC(dat = dat_c, 
                idxTrain = itrain,
                idxTest = itest,
                labelName = "professional.diagnosis", 
                featNames = featNames,
                negClassName = "FALSE", 
                posClassName = "TRUE")

## Generate the restricted permutation null
## using the random forest classifier.
##
set.seed(1234567)
resRf <- RestrictedPermAUC(dat = dat_c, 
                           idxTrain = itrain,
                           idxTest = itest,
                           nperm = nperm,
                           labelName = "professional.diagnosis", 
                           confName = "combinedConfounder",
                           featNames = featNames,
                           negClassName = "FALSE", 
                           posClassName = "TRUE",
                           verbose = FALSE,
                           parallel = TRUE)

## Get the observed AUC score
## using the logistic regression classifier.
##
set.seed(1234567)
obsGlm <- GetAUCGlm(dat = dat_c, 
                    idxTrain = itrain,
                    idxTest = itest,
                    labelName = "professional.diagnosis", 
                    featNames = featNames,
                    negClassName = "FALSE", 
                    posClassName = "TRUE")

## Generate the restricted permutation null
## using the logistic regression classifier.
##
set.seed(1234567)
resGlm <- RestrictedPermAUCGlm(dat = dat_c, 
                               idxTrain = itrain,
                               idxTest = itest,
                               nperm = nperm,
                               labelName = "professional.diagnosis", 
                               confName = "combinedConfounder",
                               featNames = featNames,
                               negClassName = "FALSE", 
                               posClassName = "TRUE",
                               verbose = FALSE,
                               parallel = TRUE)


####################################################
## Run the analyzes using the matching adjustment.
####################################################

## Get the matched training and test sets.
##
set.seed(1234567)
mtrain_dat_c <- MatchedDataByConfounder(train_dat_c, labelName = "professional.diagnosis", 
                                        confName = "combinedConfounder")
mtest_dat_c <- MatchedDataByConfounder(test_dat_c, labelName = "professional.diagnosis", 
                                       confName = "combinedConfounder")
mdat_c <- rbind(mtrain_dat_c, mtest_dat_c)
itrain <- seq(nrow(mtrain_dat_c))
itest <- max(itrain) + seq(nrow(mtest_dat_c))


set.seed(1234567)
obsRfM <- GetAUC(dat = mdat_c, 
                idxTrain = itrain,
                idxTest = itest,
                labelName = "professional.diagnosis", 
                featNames = featNames,
                negClassName = "FALSE", 
                posClassName = "TRUE")

set.seed(1234567)
resRfM <- RestrictedPermAUC(dat = mdat_c, 
                            idxTrain = itrain,
                            idxTest = itest,
                            nperm = nperm,
                            labelName = "professional.diagnosis", 
                            confName = "combinedConfounder",
                            featNames = featNames,
                            negClassName = "FALSE", 
                            posClassName = "TRUE",
                            verbose = FALSE,
                            parallel = TRUE)

set.seed(1234567)
obsGlmM <- GetAUCGlm(dat = mdat_c, 
                     idxTrain = itrain,
                     idxTest = itest,
                     labelName = "professional.diagnosis", 
                     featNames = featNames,
                     negClassName = "FALSE", 
                     posClassName = "TRUE")

set.seed(1234567)
resGlmM <- RestrictedPermAUCGlm(dat = mdat_c, 
                                idxTrain = itrain,
                                idxTest = itest,
                                nperm = nperm,
                                labelName = "professional.diagnosis", 
                                confName = "combinedConfounder",
                                featNames = featNames,
                                negClassName = "FALSE", 
                                posClassName = "TRUE",
                                verbose = FALSE,
                                parallel = TRUE)



###########################################################
## Run the analyzes using the approximate IPW adjustment.
###########################################################

## Get the augmented training and test datasets.
##
covNames <- c("age", "gender")
auxTrain <- IPWPseudoPopulation(train_dat_c, covNames)
auxTest <- IPWPseudoPopulation(test_dat_c, covNames)
itrain_dat_c <- auxTrain$rdat 
itest_dat_c <- auxTest$rdat
itrain <- seq(nrow(itrain_dat_c))
itest <- max(itrain) + seq(nrow(itest_dat_c))
idat_c <- rbind(itrain_dat_c, itest_dat_c)


set.seed(1234567)
obsRfIPWps <- GetAUC(dat = idat_c, 
                     idxTrain = itrain,
                     idxTest = itest,
                     labelName = "professional.diagnosis", 
                     featNames = featNames,
                     negClassName = "FALSE", 
                     posClassName = "TRUE")

set.seed(1234567)
resRfIPWps <- RestrictedPermAUC(dat = idat_c, 
                                idxTrain = itrain,
                                idxTest = itest,
                                nperm = nperm,
                                labelName = "professional.diagnosis", 
                                confName = "combinedConfounder",
                                featNames = featNames,
                                negClassName = "FALSE", 
                                posClassName = "TRUE",
                                verbose = FALSE,
                                parallel = TRUE)

set.seed(1234567)
obsGlmIPWps <- GetAUCGlm(dat = idat_c, 
                         idxTrain = itrain,
                         idxTest = itest,
                         labelName = "professional.diagnosis", 
                         featNames = featNames,
                         negClassName = "FALSE", 
                         posClassName = "TRUE")

set.seed(1234567)
resGlmIPWps <- RestrictedPermAUCGlm(dat = idat_c, 
                                    idxTrain = itrain,
                                    idxTest = itest,
                                    nperm = nperm,
                                    labelName = "professional.diagnosis", 
                                    confName = "combinedConfounder",
                                    featNames = featNames,
                                    negClassName = "FALSE", 
                                    posClassName = "TRUE",
                                    verbose = FALSE,
                                    parallel = TRUE)


## Save all the results.
save(obsRf, resRf, obsGlm, resGlm, 
     obsRfM, resRfM, obsGlmM, resGlmM, 
     obsRfIPWps, resRfIPWps, obsGlmIPWps, resGlmIPWps,
     file = "outputs_real_data_examples.RData", compress = TRUE)



##############################
## Plot Figures 8 and 9
##############################


hist(resRf$restrictedPermNull, probability = TRUE)
auc.u.rf <- GetUnconfoundedAUC(obsRf, resRf$restrictedPermNull)
auc.u.rf
ApproximateConfoundingPval(obsRf, resRf$restrictedPermNull)

hist(resGlm$restrictedPermNull, probability = TRUE)
auc.u.glm <- GetUnconfoundedAUC(obsGlm, resGlm$restrictedPermNull)
auc.u.glm
ApproximateConfoundingPval(obsGlm, resGlm$restrictedPermNull)


hist(resRfM$restrictedPermNull, probability = TRUE)
auc.u.rfM <- GetUnconfoundedAUC(obsRfM, resRfM$restrictedPermNull)
auc.u.rfM
ApproximateConfoundingPval(obsRfM, resRfM$restrictedPermNull)


hist(resGlmM$restrictedPermNull, probability = TRUE)
auc.u.glmM <- GetUnconfoundedAUC(obsGlmM, resGlmM$restrictedPermNull)
auc.u.glmM
ApproximateConfoundingPval(obsGlmM, resGlmM$restrictedPermNull)


hist(resRfIPWps$restrictedPermNull, probability = TRUE)
auc.u.rfIpw <- GetUnconfoundedAUC(obsRfIPWps, resRfIPWps$restrictedPermNull)
auc.u.rfIpw
ApproximateConfoundingPval(obsRfIPWps, resRfIPWps$restrictedPermNull)


hist(resGlmIPWps$restrictedPermNull, probability = TRUE)
auc.u.glmIpw <- GetUnconfoundedAUC(obsGlmIPWps, resGlmIPWps$restrictedPermNull)
auc.u.glmIpw
ApproximateConfoundingPval(obsGlmIPWps, resGlmIPWps$restrictedPermNull)


nc <- 30
xlim <- c(0.38, 0.9)
ylim <- c(0, 30)
mylwd <- 3
cl <- 1.2
cl2 <- 1.6
cm <- 1.5

## Figure 8

xaxis <- seq(0.38, 0.62, length.out = 1000)
sdensi1 <- dnorm(xaxis, 0.5, sqrt(obsGlm$approxVar["v"]))
sdensi2 <- dnorm(xaxis, 0.5, sqrt(obsGlmM$approxVar["v"]))
sdensi3 <- dnorm(xaxis, 0.5, sqrt(obsGlmIPWps$approxVar["v"]))

par(mfrow = c(3, 1), mar = c(3, 2.5, 1.5, 0.5), mgp = c(1.5, 0.5, 0))
hist(resGlm$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE,
     main = "logistic regression (no adjustment)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(a)", bty = "n", cex = cl2)
abline(v = obsGlm$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.glm, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi1, col = "red")
legend("topleft", 
       legend = c("observed AUC", "unconfounded AUC", "standard perm. null", "restricted perm. null"),
       bty = "n", text.col = c("cyan", "darkorange", "red", "blue"), cex = 1.1)
####
hist(resGlmM$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "logistic regression (matching)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(b)", bty = "n", cex = cl2)
abline(v = obsGlmM$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.glmM, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi2, col = "red")
####
hist(resGlmIPWps$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "logistic regression (approximate IPW)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(c)", bty = "n", cex = cl2)
abline(v = obsGlmIPWps$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.glmIpw, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi3, col = "red")


## Figure 9

xaxis <- seq(0.38, 0.62, length.out = 1000)
sdensi1 <- dnorm(xaxis, 0.5, sqrt(obsRf$approxVar["v"]))
sdensi2 <- dnorm(xaxis, 0.5, sqrt(obsRfM$approxVar["v"]))
sdensi3 <- dnorm(xaxis, 0.5, sqrt(obsRfIPWps$approxVar["v"]))

par(mfrow = c(3, 1), mar = c(3, 2.5, 1.5, 0.5), mgp = c(1.5, 0.5, 0))
hist(resRf$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE,
     main = "random forest (no adjustment)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(a)", bty = "n", cex = cl2)
abline(v = obsRf$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rf, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi1, col = "red")
legend("topleft", 
       legend = c("observed AUC", "unconfounded AUC", "standard perm. null", "restricted perm. null"),
       bty = "n", text.col = c("cyan", "darkorange", "red", "blue"), cex = 1.1)
####
hist(resRfM$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "random forest (matching)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(b)", bty = "n", cex = cl2)
abline(v = obsRfM$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rfM, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi2, col = "red")
####
hist(resRfIPWps$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "random forest (approximate IPW)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topright", legend = "(c)", bty = "n", cex = cl2)
abline(v = obsRfIPWps$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rfIpw, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi3, col = "red")



##############################
## Plot Figure 6
##############################

mDat$PD.status <- mDat$professional.diagnosis
mDat$PD.status[mDat$PD.status == FALSE] <- "Control"
mDat$PD.status[mDat$PD.status == TRUE] <- "PD"
mDat$gender[mDat$gender == "Male"] <- "Males"
mDat$gender[mDat$gender == "Female"] <- "Females"

cmDat <- CollapseDemo(x = mDat, subjectIdName = "healthCode")

mDatTrain <- cmDat[cmDat$datSplit == "train",]
mDatTest <- cmDat[cmDat$datSplit == "test",]


lcex <- 0.85
lcex2 <- 1.2


par(mfrow = c(2, 2), mar = c(3.5, 2.25, 2.5, 0.5) + 0.1, mgp = c(1.5, 0.5, 0))
####
hist(mDatTrain$age[mDatTrain$professional.diagnosis == FALSE], probability = FALSE,
     col = rgb(1, 0.5, 0, 0.5), xlim = c(17, 91), xlab = "age", main = "training set")
hist(mDatTrain$age[mDatTrain$professional.diagnosis == TRUE], probability = FALSE,
     col = rgb(0, 0.5, 1, 0.5), add = TRUE)
abline(v = 45, col = "red")
abline(v = 65, col = "red")
text(30, 250, "young age", cex = lcex)
text(55, 250, "middle age", cex = lcex)
text(75, 250, "senior age", cex = lcex)
legend("topright", legend = c("Control", "PD"), 
       text.col = c(rgb(1, 0.5, 0, 1), rgb(0, 0.5, 1, 1)), bty = "n")
mtext("(a)", side = 3, at = 18, cex = lcex2)
####
hist(mDatTest$age[mDatTest$professional.diagnosis == FALSE], probability = FALSE,
     col = rgb(1, 0.5, 0, 0.5), xlim = c(17, 91), xlab = "age", main = "test set")
hist(mDatTest$age[mDatTest$professional.diagnosis == TRUE], probability = FALSE,
     col = rgb(0, 0.5, 1, 0.5), add = TRUE)
abline(v = 45, col = "red")
abline(v = 65, col = "red")
text(30, 140, "young age", cex = lcex)
text(55, 140, "middle age", cex = lcex)
text(75, 140, "senior age", cex = lcex)
legend("topright", legend = c("Control", "PD"), 
       text.col = c(rgb(1, 0.5, 0, 1), rgb(0, 0.5, 1, 1)), bty = "n")
mtext("(b)", side = 3, at = 18, cex = lcex2)
####
par(mar = c(2.5, 2.25, 1.5, 0.5) + 0.1)
####
mosaicplot(PD.status ~ gender, 
           data = mDatTrain, 
           shade = FALSE, 
           color = TRUE,
           main = "training set",
           xlab = "DISEASE STATUS",
           ylab = "GENDER",
           cex.axis = 0.9)
mtext("(c)", side = 3, at = 0.05, cex = lcex2)
####
mosaicplot(PD.status ~ gender, 
           data = mDatTest, 
           shade = FALSE, 
           color = TRUE,
           main = "test set",
           xlab = "DISEASE STATUS",
           ylab = "GENDER",
           cex.axis = 0.9)
mtext("(d)", side = 3, at = 0.05, cex = lcex2)
par(mar = c(5, 4, 4, 2) + 0.1)



#############################
## Plot Figure 10
#############################

train_dat_c$professional.diagnosis2 <- train_dat_c$professional.diagnosis
train_dat_c$professional.diagnosis2[train_dat_c$professional.diagnosis == "FALSE"] <- "Control"
train_dat_c$professional.diagnosis2[train_dat_c$professional.diagnosis == "TRUE"] <- "PD"
test_dat_c$professional.diagnosis2 <- test_dat_c$professional.diagnosis
test_dat_c$professional.diagnosis2[test_dat_c$professional.diagnosis == "FALSE"] <- "Control"
test_dat_c$professional.diagnosis2[test_dat_c$professional.diagnosis == "TRUE"] <- "PD"

itrain_dat_c$professional.diagnosis2 <- itrain_dat_c$professional.diagnosis
itrain_dat_c$professional.diagnosis2[itrain_dat_c$professional.diagnosis == "FALSE"] <- "Control"
itrain_dat_c$professional.diagnosis2[itrain_dat_c$professional.diagnosis == "TRUE"] <- "PD"
itest_dat_c$professional.diagnosis2 <- itest_dat_c$professional.diagnosis
itest_dat_c$professional.diagnosis2[itest_dat_c$professional.diagnosis == "FALSE"] <- "Control"
itest_dat_c$professional.diagnosis2[itest_dat_c$professional.diagnosis == "TRUE"] <- "PD"

mtrain_dat_c$professional.diagnosis2 <- mtrain_dat_c$professional.diagnosis
mtrain_dat_c$professional.diagnosis2[mtrain_dat_c$professional.diagnosis == "FALSE"] <- "Control"
mtrain_dat_c$professional.diagnosis2[mtrain_dat_c$professional.diagnosis == "TRUE"] <- "PD"
mtest_dat_c$professional.diagnosis2 <- mtest_dat_c$professional.diagnosis
mtest_dat_c$professional.diagnosis2[mtest_dat_c$professional.diagnosis == "FALSE"] <- "Control"
mtest_dat_c$professional.diagnosis2[mtest_dat_c$professional.diagnosis == "TRUE"] <- "PD"


par(mfrow = c(3, 2), mar = c(3, 1, 3, 0.1), mgp = c(1.75, 0.75, 0))
mosaicplot(train_dat_c$professional.diagnosis2 ~ train_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "no adjustment, training set")
mtext(side = 3, "(a)", at = 0, line = -1)
mosaicplot(test_dat_c$professional.diagnosis2 ~ test_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "no adjustment, test set")
mtext(side = 3, "(b)", at = 0, line = -1)
####
mosaicplot(itrain_dat_c$professional.diagnosis2 ~ itrain_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "approx. IPW, training set")
mtext(side = 3, "(c)", at = 0, line = -1)
mosaicplot(itest_dat_c$professional.diagnosis2 ~ itest_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "approx. IPW, test set")
mtext(side = 3, "(d)", at = 0, line = -1)
####
mosaicplot(mtrain_dat_c$professional.diagnosis ~ mtrain_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "matching, training set")
mtext(side = 3, "(e)", at = 0, line = -1)
mosaicplot(mtest_dat_c$professional.diagnosis2 ~ mtest_dat_c$combinedConfounder, las = 2,
           ylab = "combined confounder", xlab = "disease status", main = "matching, test set")
mtext(side = 3, "(f)", at = 0, line = -1)



##################################################
## Run the random forest analysis for a feature
## set generated by a deep learning model
## (to be used in Supplementary Figure S2).
##################################################

require(synapseClient)
synapseLogin("echaibub", "2011synapse")

load(getFileLocation(synGet("syn12182549")))
feats <- read.csv(getFileLocation(synGet("syn10949406")), header = TRUE)

dat <- data.frame(mDat, feats[match(mDat$recordId, feats$recordId),])
dat$age <- as.numeric(dat$age)
dat <- na.omit(dat)
idx_train <- which(dat$datSplit == "train")
idx_test <- which(dat$datSplit == "test")
dat_train <- dat[idx_train,]
dat_test <- dat[idx_test,]
featNames <- colnames(feats)[-1]

cat("collapse features", "\n")
train_dat_c <- CollapseFeaturesMedian(x = dat_train, 
                                      labelName = "professional.diagnosis", 
                                      covNames = c("age", "gender"), 
                                      subjectIdName = "healthCode", 
                                      featNames = featNames)
test_dat_c <- CollapseFeaturesMedian(x = dat_test, 
                                     labelName = "professional.diagnosis", 
                                     covNames = c("age", "gender"), 
                                     subjectIdName = "healthCode", 
                                     featNames = featNames)


train_dat_c$dage <- GetDiscretizedAge(train_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                      levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))
test_dat_c$dage <- GetDiscretizedAge(test_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                     levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))

train_dat_c$combinedConfounder <- paste(train_dat_c$dage, train_dat_c$gender, sep = "_")
test_dat_c$combinedConfounder <- paste(test_dat_c$dage, test_dat_c$gender, sep = "_")

dat_c <- rbind(train_dat_c, test_dat_c)

nperm <- 1e+4

########################
## no adjustment 
########################

itrain <- seq(nrow(train_dat_c))
itest <- max(itrain) + seq(nrow(test_dat_c))

set.seed(1234567)
obsRf <- GetAUC(dat = dat_c, 
                idxTrain = itrain,
                idxTest = itest,
                labelName = "professional.diagnosis", 
                featNames = featNames,
                negClassName = "FALSE", 
                posClassName = "TRUE")

set.seed(1234567)
resRf <- RestrictedPermAUC(dat = dat_c, 
                           idxTrain = itrain,
                           idxTest = itest,
                           nperm = nperm,
                           labelName = "professional.diagnosis", 
                           confName = "combinedConfounder",
                           featNames = featNames,
                           negClassName = "FALSE", 
                           posClassName = "TRUE",
                           verbose = FALSE,
                           parallel = TRUE)


##########################
## Matching adjustment
##########################

set.seed(1234567)
mtrain_dat_c <- MatchedDataByConfounder(train_dat_c, labelName = "professional.diagnosis", 
                                        confName = "combinedConfounder")
mtest_dat_c <- MatchedDataByConfounder(test_dat_c, labelName = "professional.diagnosis", 
                                       confName = "combinedConfounder")
mdat_c <- rbind(mtrain_dat_c, mtest_dat_c)
itrain <- seq(nrow(mtrain_dat_c))
itest <- max(itrain) + seq(nrow(mtest_dat_c))


set.seed(1234567)
obsRfM <- GetAUC(dat = mdat_c, 
                 idxTrain = itrain,
                 idxTest = itest,
                 labelName = "professional.diagnosis", 
                 featNames = featNames,
                 negClassName = "FALSE", 
                 posClassName = "TRUE")

set.seed(1234567)
resRfM <- RestrictedPermAUC(dat = mdat_c, 
                            idxTrain = itrain,
                            idxTest = itest,
                            nperm = nperm,
                            labelName = "professional.diagnosis", 
                            confName = "combinedConfounder",
                            featNames = featNames,
                            negClassName = "FALSE", 
                            posClassName = "TRUE",
                            verbose = FALSE,
                            parallel = TRUE)


################################
## Approximate IPW adjustment
################################

covNames <- c("age", "gender")
auxTrain <- IPWPseudoPopulation(train_dat_c, covNames)
auxTest <- IPWPseudoPopulation(test_dat_c, covNames)

itrain_dat_c <- auxTrain$rdat 
itest_dat_c <- auxTest$rdat
itrain <- seq(nrow(itrain_dat_c))
itest <- max(itrain) + seq(nrow(itest_dat_c))
idat_c <- rbind(itrain_dat_c, itest_dat_c)


set.seed(1234567)
obsRfIPWps <- GetAUC(dat = idat_c, 
                     idxTrain = itrain,
                     idxTest = itest,
                     labelName = "professional.diagnosis", 
                     featNames = featNames,
                     negClassName = "FALSE", 
                     posClassName = "TRUE")


set.seed(1234567)
resRfIPWps <- RestrictedPermAUC(dat = idat_c, 
                                idxTrain = itrain,
                                idxTest = itest,
                                nperm = nperm,
                                labelName = "professional.diagnosis", 
                                confName = "combinedConfounder",
                                featNames = featNames,
                                negClassName = "FALSE", 
                                posClassName = "TRUE",
                                verbose = FALSE,
                                parallel = TRUE)


save(obsRf, resRf, obsRfM, resRfM, obsRfIPWps, resRfIPWps, 
     file = "outputs_real_data_example_deep_model_feats.RData", compress = TRUE)


######################################
## Plot Supplementary Figure S2
######################################

auc.u.rf <- GetUnconfoundedAUC(obsRf, resRf$restrictedPermNull)
auc.u.rfM <- GetUnconfoundedAUC(obsRfM, resRfM$restrictedPermNull)
auc.u.rfIpw <- GetUnconfoundedAUC(obsRfIPWps, resRfIPWps$restrictedPermNull)

nc <- 30
xlim <- c(0.38, 0.9)
ylim <- c(0, 40)
mylwd <- 3
cl <- 1.2
cl2 <- 1.6
cm <- 1.5
ca <- 1.2

xaxis <- seq(0.38, 0.62, length.out = 1000)
sdensi1 <- dnorm(xaxis, 0.5, sqrt(obsRf$approxVar["v"]))
sdensi2 <- dnorm(xaxis, 0.5, sqrt(obsRfM$approxVar["v"]))
sdensi3 <- dnorm(xaxis, 0.5, sqrt(obsRfIPWps$approxVar["v"]))

par(mfrow = c(3, 1), mar = c(3, 2.5, 1.5, 0.5), mgp = c(1.5, 0.5, 0))
hist(resRf$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE,
     main = "random forest (no adjustment)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topleft", legend = "(a)", bty = "n", cex = cl2)
abline(v = obsRf$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rf, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi1, col = "red")
####
hist(resRfM$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "random forest (matching)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topleft", legend = "(b)", bty = "n", cex = cl2)
abline(v = obsRfM$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rfM, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi2, col = "red")
####
hist(resRfIPWps$restrictedPermNull, xlim = xlim, ylim = ylim,  xlab = "AUC", probability = TRUE, 
     main = "random forest (approximate IPW)", col = rgb(0, 0, 1, 1), nclass = nc, cex.main = cm,
     cex.axis = ca, cex.lab = cl)
legend("topleft", legend = "(c)", bty = "n", cex = cl2)
abline(v = obsRfIPWps$aucObs, col = "cyan", lwd = mylwd)
abline(v = auc.u.rfIpw, col = "darkorange", lwd = mylwd)
lines(xaxis, sdensi3, col = "red")



