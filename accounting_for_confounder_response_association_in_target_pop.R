
source("utility_functions_kdd.R")

## Simulate the development dataset.

n <- 10000
idxTrain <- seq(1, n/2)
idxTest <- seq(n/2 + 1, n)
nfeat <- 10
featNames <- paste("f", seq(nfeat), sep = "")

set.seed(1234567)
dat1 <- GenerateData(n,
                     nfeat = nfeat, 
                     p11 = 0.4,
                     p10 = 0.1,
                     p01 = 0.1,
                     p00 = 0.4,
                     alpha = c(0.25, 1),
                     rho = 0.5,
                     binVarNames = c("disease", "gender"))
dat1$disease <- factor(dat1$disease, labels = c("control", "case"))
dat1$gender <- factor(dat1$gender, labels = c("female", "male"))

## Gender by disease label table
## from the developement data.
##
table(dat1$gender, dat1$disease)

## Frequency of the gender by disease label
## classes in the subsets from the development
## data that will form the baseline training
## and test sets.
##
## Note that the number of control males (986) 
## is the "limiting factor" here. 
##
## Because in the target population we have that
## half of the controls are male, we need half of
## the controls to be males in the baseline
## training and test sets. Hence, we assign 489
## control males to baseline training set and 489 
## control males to the baseline test set. 
##
## Now, since the controls make 66.66% of the data
## and we want 11.11% to be female cases and 22.22% 
## to be male cases, we need to sample 163 and 326
## female cases and male cases, respectively.
##
freqs <- matrix(c(489, 163,
                  489, 326), 2, 2, byrow = TRUE)
colnames(freqs) <- c("control", "case")
rownames(freqs) <- c("female", "male")

## Subsample the baseline training and test sets from the
## development data (according to the frequencies that 
## were set up above).
##
auxBase <- GetBaselineTrainTestSets(dat1, freqs, respName, confName)
baseTrain <- auxBase$datTrain
baseTest <- auxBase$datTest

baseDat <- rbind(baseTrain, baseTest)
idxTrain <- seq(nrow(baseTrain))
idxTest <- max(idxTrain) + seq(nrow(baseTest))

nperm <- 1e+4

## Run the restricted permutations on the baseline
## training and test sets.
##
set.seed(12345)
resB <- RestrictedPermAUC(baseDat, 
                          idxTrain,
                          idxTest,
                          nperm = nperm,
                          labelName = "disease", 
                          confName = "gender",
                          featNames = featNames,
                          negClassName = "control", 
                          posClassName = "case",
                          verbose = FALSE,
                          parallel = TRUE)


## Split the development data into the 
## development training set and 
## development test set.
##
ntest <- sum(freqs)
aux <- TrainTestStratifiedSplit(dat1, ntest = ntest, respName = "disease", confName = "gender")
aux$trainNumbers
aux$testNumbers
table(dat1[aux$idxTrain, "gender"], dat1[aux$idxTrain, "disease"])
table(dat1[aux$idxTest, "gender"], dat1[aux$idxTest, "disease"])

## Calculate the observed AUC score (and the 
## approximated variance for the standard
## permutation null distribution).
##
set.seed(12345)
obs <- GetAUC(dat1,
              idxTrain = aux$idxTrain,
              idxTest = aux$idxTest,
              labelName = "disease", 
              featNames = featNames,
              negClassName = "control", 
              posClassName = "case")

## Generate the restricted permutation null distribution
## from the development data.
##
set.seed(12345)
resD <- RestrictedPermAUC(dat1, 
                          idxTrain = aux$idxTrain,
                          idxTest = aux$idxTest,
                          nperm = nperm,
                          labelName = "disease", 
                          confName = "gender",
                          featNames = featNames,
                          negClassName = "control", 
                          posClassName = "case",
                          verbose = FALSE,
                          parallel = TRUE)

## Save the results.
##
save(resB, resD, file = "output_restricted_perms_Baseline_Development_sets.RData", compress = TRUE)



## Generate a matrix representing the a priori
## knowledge about the association between disease
## and gender in the target population (1/3 of the 
## population has the disease, and the disease is 
## two times more common in men than in women).
##
## (This matrix is used to generate panel a in Figure 11.)
##
popInt <- matrix(c("control", "female",
                   "control", "female",
                   "control", "female",
                   "control", "male",
                   "control", "male",
                   "control", "male",
                   "case", "female",
                   "case", "male",
                   "case", "male"), 9, 2, byrow = TRUE)
popInt <- data.frame(popInt)
names(popInt) <- c("disease", "gender")
popInt$disease <- factor(popInt$disease, levels = c("control", "case"))



#########################
## Plot panels a and b
#########################

mycex <- 1.2
ca <- 1.3
cm <- 1.7
ca2 <- 1.5

par(mfrow = c(1, 2))
mosaicplot(disease ~ gender, data = popInt, cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "population of interest", cex = cm, line = 1.5)
mtext(side = 3, "(a)", line = -0.5, at = 0, cex = cm)
text(0.35, 0.7, c("(33.33%)"), cex = mycex)
text(0.35, 0.2, c("(33.33%)"), cex = mycex)
text(0.83, 0.8, c("(11.11%)"), cex = mycex)
text(0.83, 0.3, c("(22.22%)"), cex = mycex)
####
mosaicplot(disease ~ gender, data = dat1, cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "development population", cex = cm, line = 1.5)
mtext(side = 3, "(b)", line = -0.5, at = 0, cex = cm)
text(0.3, 0.65, c("4082"), cex = mycex)
text(0.3, 0.575, c("(40.82%)"), cex = mycex)
text(0.3, 0.1, c("986"), cex = mycex)
text(0.3, 0.025, c("(9.86%)"), cex = mycex)
text(0.75, 0.9, c("978"), cex = mycex)
text(0.75, 0.825, c("(9.78%)"), cex = mycex)
text(0.75, 0.4, c("3954"), cex = mycex)
text(0.75, 0.325, c("(39.54%)"), cex = mycex)
par(mfrow = c(1, 1))


#########################
## Plot panels c and d
#########################

## Get the numbers for baseline training and test sets.
##
table(baseTrain$gender, baseTrain$disease)
round(100*table(baseTrain$gender, baseTrain$disease)/nrow(baseTrain), 2)
table(baseTest$gender, baseTest$disease)
round(100*table(baseTest$gender, baseTest$disease)/nrow(baseTest), 2)

par(mfrow = c(1, 2))
mosaicplot(disease ~ gender, data = baseTrain, cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "baseline train set", cex = cm, line = 1.5)
mtext(side = 3, "(c)", line = -0.5, at = 0, cex = cm)
text(0.35, 0.74, c("489"), cex = mycex)
text(0.35, 0.66, c("(33.33%)"), cex = mycex)
text(0.35, 0.24, c("489"), cex = mycex)
text(0.35, 0.16, c("(33.33%)"), cex = mycex)
text(0.83, 0.84, c("163"), cex = mycex)
text(0.83, 0.76, c("(11.11%)"), cex = mycex)
text(0.83, 0.34, c("326"), cex = mycex)
text(0.83, 0.26, c("(22.22%)"), cex = mycex)
####
mosaicplot(disease ~ gender, data = baseTest, cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "baseline test set", cex = cm, line = 1.5)
mtext(side = 3, "(d)", line = -0.5, at = 0, cex = cm)
text(0.35, 0.74, c("489"), cex = mycex)
text(0.35, 0.66, c("(33.33%)"), cex = mycex)
text(0.35, 0.24, c("489"), cex = mycex)
text(0.35, 0.16, c("(33.33%)"), cex = mycex)
text(0.83, 0.84, c("163"), cex = mycex)
text(0.83, 0.76, c("(11.11%)"), cex = mycex)
text(0.83, 0.34, c("326"), cex = mycex)
text(0.83, 0.26, c("(22.22%)"), cex = mycex)
par(mfrow = c(1, 1))



#########################
## Plot panels f and g
#########################

## Get the numbers for development training and test sets.
##
table(dat1[aux$idxTrain, "gender"], dat1[aux$idxTrain, "disease"])
round(100*table(dat1[aux$idxTrain, "gender"], dat1[aux$idxTrain, "disease"])/length(aux$idxTrain), 2)
table(dat1[aux$idxTest, "gender"], dat1[aux$idxTest, "disease"])
round(100*table(dat1[aux$idxTest, "gender"], dat1[aux$idxTest, "disease"])/length(aux$idxTest), 2)

par(mfrow = c(1, 2))
mosaicplot(disease ~ gender, data = dat1[aux$idxTrain,], cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "development train set", cex = cm, line = 1.5)
mtext(side = 3, "(f)", line = -0.5, at = 0, cex = cm)
text(0.3, 0.65, c("3483"), cex = mycex)
text(0.3, 0.575, c("(40.82%)"), cex = mycex)
text(0.3, 0.1, c("841"), cex = mycex)
text(0.3, 0.025, c("(9.86%)"), cex = mycex)
text(0.75, 0.9, c("835"), cex = mycex)
text(0.75, 0.825, c("(9.79%)"), cex = mycex)
text(0.75, 0.4, c("3374"), cex = mycex)
text(0.75, 0.325, c("(39.54%)"), cex = mycex)
####
mosaicplot(disease ~ gender, data = dat1[aux$idxTest,], cex.axis = ca, 
           main = "", xlab = "", ylab = "")
mtext(side = 1, "disease", cex = ca2, line = 1.5)
mtext(side = 2, "gender", cex = ca2, line = 1.5)
mtext(side = 3, "development test set", cex = cm, line = 1.5)
mtext(side = 3, "(g)", line = -0.5, at = 0, cex = cm)
text(0.3, 0.65, c("599"), cex = mycex)
text(0.3, 0.575, c("(40.83%)"), cex = mycex)
text(0.3, 0.1, c("145"), cex = mycex)
text(0.3, 0.025, c("(9.88%)"), cex = mycex)
text(0.75, 0.9, c("143"), cex = mycex)
text(0.75, 0.825, c("(9.75%)"), cex = mycex)
text(0.75, 0.4, c("580"), cex = mycex)
text(0.75, 0.325, c("(39.54%)"), cex = mycex)
par(mfrow = c(1, 1))


####################
## Plot panel e
####################

load("output_restricted_perms_Baseline_Development_sets.RData")

## Compute the unconfounded AUC estimate (relative 
## to the target distribution).
##
auc.0 <- obs$aucObs
a.D <- mean(resD$restrictedPermNull)
s.D <- sd(resD$restrictedPermNull)
a.B <- mean(resB$restrictedPermNull)
s.B <- sd(resB$restrictedPermNull)
auc.u <- (auc.0 - a.D) * (s.B/s.D) + a.B
auc.u

## Compute the unconfounded AUC estimate (relative 
## to the standard permutation null distribution).
##
avarStandNull <- as.numeric(obs$approxVar["v"])
s.2.stars <- sqrt(avarStandNull)
auc.u.2 <- (auc.0 - a.D) * (s.2.stars/s.D) + 0.5
auc.u.2

nc <- 30
myxlim <- c(0.4, 0.9)
myylim <- NULL
mylwd <- 2
a <- 1

xaxis <- seq(0.43, 0.57, length.out = 1000)
densi <- dnorm(xaxis, 0.5, sqrt(obs$approxVar["v"]))

hist(resD$restrictedPermNull, probability = TRUE, col = rgb(0, 0, 1, 0.25), nclass = nc, xlim = myxlim, 
     xlab = "AUC", main = "", ylim = myylim)
hist(resB$restrictedPermNull, probability = TRUE, col = rgb(0, 1, 0, 0.25), nclass = nc, xlim = myxlim, 
     xlab = "AUC", main = "", ylim = myylim, add = TRUE)
abline(v = auc.0, col = "cyan", lwd = mylwd)
abline(v = auc.u.2, col = "darkorange", lwd = mylwd)
abline(v = auc.u, col = "darkgreen", lwd = mylwd)
lines(xaxis, densi, col = "red")





