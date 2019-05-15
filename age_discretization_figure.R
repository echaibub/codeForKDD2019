
## Source utility functions.
##
source("utility_functions_kdd.R")

##################################
## Additional utility functions.
##################################

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


########################################################
########################################################
########################################################

#############################
## Load and shape the data.
#############################

require(synapseClient)
synapseLogin()

nperm <- 1e+4

## load metadata
load(getFileLocation(synGet("syn12182549")))

## load features
feats <- read.csv(getFileLocation(synGet("syn11027900")), header = TRUE)

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

train_dat_c$dage2 <- GetDiscretizedAge(train_dat_c, breaks = c(17, 58, 99), nlevels = 2, 
                                       levelNames = c("YoungAge", "SeniorAge"))
test_dat_c$dage2 <- GetDiscretizedAge(test_dat_c, breaks = c(17, 58, 99), nlevels = 2, 
                                      levelNames = c("YoungAge", "SeniorAge"))

train_dat_c$dage3 <- GetDiscretizedAge(train_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                       levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))
test_dat_c$dage3 <- GetDiscretizedAge(test_dat_c, breaks = c(17, 44, 65, 99), nlevels = 3, 
                                      levelNames = c("YoungAge", "MiddleAge", "SeniorAge"))

train_dat_c$dage4 <- GetDiscretizedAge(train_dat_c, breaks = c(17, 35, 50, 65, 99), nlevels = 4, 
                                       levelNames = paste("a", 1:4, sep = ""))
test_dat_c$dage4 <- GetDiscretizedAge(test_dat_c, breaks = c(17, 35, 50, 65, 99), nlevels = 4, 
                                      levelNames = paste("a", 1:4, sep = ""))

train_dat_c$dage5 <- GetDiscretizedAge(train_dat_c, breaks = c(17, 30, 45, 60, 75, 99), nlevels = 5, 
                                       levelNames = paste("a", 1:5, sep = ""))
test_dat_c$dage5 <- GetDiscretizedAge(test_dat_c, breaks = c(17, 30, 45, 60, 75, 99), nlevels = 5, 
                                      levelNames = paste("a", 1:5, sep = ""))

train_dat_c$combinedConfounder <- paste(train_dat_c$dage, train_dat_c$gender, sep = "_")
test_dat_c$combinedConfounder <- paste(test_dat_c$dage, test_dat_c$gender, sep = "_")

dat_c <- rbind(train_dat_c, test_dat_c)
itrain <- seq(nrow(train_dat_c))
itest <- seq(nrow(test_dat_c)) + nrow(train_dat_c)


#########################################
## Generate the restricted permutation
## null distributions for distinct
## age discretizations.
#########################################

myseed <- 1234567890

cat("res2", "\n")
set.seed(myseed)
res2 <- RestrictedPermAUC(dat = dat_c, 
                          idxTrain = itrain,
                          idxTest = itest,
                          nperm = nperm,
                          labelName = "professional.diagnosis", 
                          confName = "dage2",
                          featNames = featNames,
                          negClassName = "FALSE", 
                          posClassName = "TRUE",
                          verbose = FALSE,
                          parallel = TRUE)

cat("res3", "\n")
set.seed(myseed)
res3 <- RestrictedPermAUC(dat = dat_c, 
                          idxTrain = itrain,
                          idxTest = itest,
                          nperm = nperm,
                          labelName = "professional.diagnosis", 
                          confName = "dage3",
                          featNames = featNames,
                          negClassName = "FALSE", 
                          posClassName = "TRUE",
                          verbose = FALSE,
                          parallel = TRUE)

cat("res4", "\n")
set.seed(myseed)
res4 <- RestrictedPermAUC(dat = dat_c, 
                          idxTrain = itrain,
                          idxTest = itest,
                          nperm = nperm,
                          labelName = "professional.diagnosis", 
                          confName = "dage4",
                          featNames = featNames,
                          negClassName = "FALSE", 
                          posClassName = "TRUE",
                          verbose = FALSE,
                          parallel = TRUE) 

cat("res5", "\n")
set.seed(myseed)
res5 <- RestrictedPermAUC(dat = dat_c, 
                          idxTrain = itrain,
                          idxTest = itest,
                          nperm = nperm,
                          labelName = "professional.diagnosis", 
                          confName = "dage5",
                          featNames = featNames,
                          negClassName = "FALSE", 
                          posClassName = "TRUE",
                          verbose = FALSE,
                          parallel = TRUE)



#########################
## Plot Figure 7
#########################

xlim <- c(0.5, 0.76)
cl <- 1.2
cl2 <- 1.6
nc <- 30
mylwd <- 2


par(mfrow = c(2, 2), mar = c(3.5, 3, 1, 1), mgp = c(1.5, 0.5, 0))
hist(res2$restrictedPermNull, xlim = xlim, xlab = "AUC", main = "2 levels",
     col = rgb(0, 0, 1, 0.5), nclass = nc)
abline(v = mean(res2$restrictedPermNull), col = "red", lwd = mylwd)
legend("topleft", legend = "(a)", bty = "n", cex = cl2)
####
hist(res3$restrictedPermNull, xlim = xlim, xlab = "AUC", main = "3 levels",
     col = rgb(0, 0, 1, 0.5), nclass = nc)
abline(v = mean(res3$restrictedPermNull), col = "red", lwd = mylwd)
legend("topleft", legend = "(b)", bty = "n", cex = cl2)
####
hist(res4$restrictedPermNull, xlim = xlim, xlab = "AUC", main = "4 levels",
     col = rgb(0, 0, 1, 0.5), nclass = nc)
abline(v = mean(res4$restrictedPermNull), col = "red", lwd = mylwd)
legend("topleft", legend = "(c)", bty = "n", cex = cl2)
####
hist(res5$restrictedPermNull, xlim = xlim, xlab = "AUC", main = "5 levels",
     col = rgb(0, 0, 1, 0.5), nclass = nc)
abline(v = mean(res5$restrictedPermNull), col = "red", lwd = mylwd)
legend("topleft", legend = "(d)", bty = "n", cex = cl2)

