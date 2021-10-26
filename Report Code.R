####### ORIGINAL CODE - BEFORE EXTENSION #######

## The code below uses the following datasets: salta_data.tab, script1_recording.R, script2_recording.R 
## that can be found here: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/24896

########################################################################
# File name: "matching_script.R"
# Goal: Estimate effect of e-voting using matching
# Dependency: "datamatch.Rdata" 
########################################################################

## Loading necessary libraries 
library(MatchIt)
library(rbounds)
#Loading the data
load("datamatch.Rdata")
#Specifying outcomes in the dataset
outcomes <- datamatch[10:18]
#Getting the names for the outcome variables
outcomes.lbls <- names(outcomes)
#Specifying the number of outcome variables
n.outcomes <- dim(outcomes)[2]

#___________________________________________________________________________#

# Dropping observations with missing values in covariates

datamatch[, 10:18][is.na(datamatch[, 10:18]) == "TRUE"] <- 99999
datamatch <- na.omit(datamatch)

#__________________________ Table 2, pre-matching balance __________________________#

#Defining the units that did e-voting (treatment units)
EV <- datamatch[2]
#Defining the covariates we want to match on
covariates <- datamatch[c("age.group", "educ", "white.collar", "not.full.time", "male", "tech", "pol.info")]
covariate.lbls <- names(covariates)
#Counting the number of covariates
n.covariates <- dim(covariates)[2]
#Creating a table with the proportions of e-voting and traditional voting, their difference, and p-values
tab2.pre <- matrix(NA, nrow = n.covariates, ncol = 4)
rownames(tab2.pre) <- covariate.lbls
colnames(tab2.pre) <- c("ev", "tv", "diff", "pvalue")

tab2.pre[, 1:2] <- cbind(apply(covariates[EV == 1,], 2, mean), apply(covariates[EV == 0,], 2, mean))
tab2.pre[, 3] <- tab2.pre[, 1] - tab2.pre[, 2]

for (i in c(1, 2, 6, 7)){
  tab2.pre[i, 4] <- ks.boot(covariates[, i][EV == 1], covariates[, i][EV == 0], nboots = 500)$ks.boot.pvalue
}
for (i in c(3, 4, 5)){
  tab2.pre[i, 4] <- prop.test(table(covariates[, i], EV$EV), n = apply(table(covariates[,i],EV$EV),2, sum))$p.value
}

#__________________________ Matching (with MatchIt) ________________________#

print("Matching")

set.seed(36466)

m.out <- matchit(EV ~ age.group + I(age.group^2) + I(age.group^3) + age.group:educ + age.group:tech + educ + I(educ^2) + tech + I(tech^2) + pol.info + educ:pol.info + age.group:pol.info + tech:pol.info + white.collar + not.full.time + male, caliper = 0.05, data = datamatch, method = "nearest", verbose = "TRUE")

print("Balance Improvement")
#summary of balance
print(summary(m.out))

#___________________________________________________________________________#

# matched sample

datamatched <- match.data(m.out)
datamatched[datamatched == 99999] <- NA

save(datamatched, file = "datamatched.Rdata")

#__________________________ Table 2, post-matching _________________________#

EV.post <- datamatched[2]

covariates.post <- datamatched[, covariate.lbls]

tab2.post <- matrix(NA, nrow = n.covariates, ncol = 4)
rownames(tab2.post) <- covariate.lbls
colnames(tab2.post) <- c("ev", "tv", "diff", "pvalue")

tab2.post[, 1:2] <- cbind(apply(covariates.post[EV.post == 1, ], 2, mean), apply(covariates.post[EV.post == 0,], 2, mean))
tab2.post[, 3] <- tab2.post[, 1] - tab2.post[, 2]
for (i in c(1, 2, 6 , 7)){
  tab2.post[i, 4]<-ks.boot(covariates.post[,i][EV.post==1],covariates.post[,i][EV.post==0], nboots = 500)$ks.boot.pvalue
}
for (i in c(3, 4, 5)){
  tab2.post[i, 4] <- prop.test(table(covariates.post[, i], EV.post$EV), n = apply(table(covariates.post[, i], EV.post$EV),2 , sum))$p.value
}

tab2 <- cbind(tab2.pre, tab2.post)
tab2[3:5, c(1:3, 5:7)] <- tab2[3:5, c(1:3, 5:7)] * 100

### Table 2 ###

print(tab2, digits = 2)
#__________________________ Table 3, pre-matching __________________________#

datamatch[datamatch == 99999] <- NA

outcomes.pre <- datamatch[10:18]

tab3.pre <- matrix(NA,nrow = n.outcomes,ncol = 5)
rownames(tab3.pre) <- outcomes.lbls
colnames(tab3.pre) <- c("N", "prop.ev", "prop.tv", "diff", "pvalue")

for (i in 1:n.outcomes) {
  tab3.pre[i, 1] <- length(na.omit(outcomes.pre[, i]))
  tab3.pre[i, 2:3] <- rev(prop.table(table(outcomes.pre[,i],datamatch$EV),2)[2,])*100
  tab3.pre[i, 4] <- tab3.pre[i, 2] - tab3.pre[i, 3]	
  tab3.pre[i, 5] <- prop.test(table(outcomes.pre[, i], datamatch$EV)[2, ], n = apply(table(outcomes.pre[, i], datamatch$EV), 2, sum))$p.value
}

datamatch[, 10:18][is.na(datamatch[, 10:18]) == "TRUE"] <- 99999


#__________________________ Table 3, post-matching _________________________#

outcomes.post <- datamatched[10:18]

tab3.post <- matrix(NA, nrow = n.outcomes, ncol = 5)
rownames(tab3.post) <- outcomes.lbls
colnames(tab3.post) <- c("N", "prop.ev", "prop.tv", "diff", "pvalue")

for (i in 1:n.outcomes) {
  tab3.post[i, 1] <- length(na.omit(outcomes.post[, i]))
  tab3.post[i, 2:3] <- rev(prop.table(table(outcomes.post[, i], datamatched$EV), 2)[2, ]) * 100
  tab3.post[i, 4] <- tab3.post[i, 2] - tab3.post[i, 3]	
  tab3.post[i, 5] <- prop.test(table(outcomes.post[, i], datamatched$EV)[2, ], n = apply(table(outcomes.post[, i], datamatched$EV), 2, sum))$p.value
}

tab3 <- cbind(tab3.pre, tab3.post)

tab3 <- tab3[rev(order(tab3[, 9])), ]

### Table 3 ###

print(tab3, digits = 4)

####### EXTENSION CODE - GENETIC MATCHING #######

#__________________________ Table 2, pre-matching balance __________________________#

#Defining the units that did e-voting (treatment units)
EV <- datamatch[2]
#Defining the covariates we want to match on
covariates <- datamatch[c("age.group", "educ", "white.collar", "not.full.time", "male", "tech", "pol.info")]
covariate.lbls <- names(covariates)
#Counting the number of covariates
n.covariates <- dim(covariates)[2]
#Creating a table with the proportions of e-voting and traditional voting, their difference, and p-values
tab2.pre <- matrix(NA, nrow = n.covariates, ncol = 4)
rownames(tab2.pre) <- covariate.lbls
colnames(tab2.pre) <- c("ev", "tv", "diff", "pvalue")

tab2.pre[, 1:2] <- cbind(apply(covariates[EV == 1,], 2, mean), apply(covariates[EV == 0,], 2, mean))
tab2.pre[, 3] <- tab2.pre[, 1] - tab2.pre[, 2]

for (i in c(1, 2, 6, 7)){
  tab2.pre[i, 4] <- ks.boot(covariates[, i][EV == 1], covariates[, i][EV == 0], nboots = 500)$ks.boot.pvalue
}
for (i in c(3, 4, 5)){
  tab2.pre[i, 4] <- prop.test(table(covariates[, i], EV$EV), n = apply(table(covariates[,i],EV$EV),2, sum))$p.value
}

#__________________________ Genetic Matching (with MatchIt) ________________________#

print("Matching")

set.seed(36577)

gm.out <- matchit(EV ~ age.group + I(age.group^2) + I(age.group^3) + age.group:educ + age.group:tech + educ + I(educ^2) + tech + I(tech^2) + pol.info + educ:pol.info + age.group:pol.info + tech:pol.info + white.collar + not.full.time + male, data = datamatch, method = "genetic", verbose = "TRUE")

print("Balance Improvement")
print(summary(gm.out))

#___________________________________________________________________________#

# matched sample

datamatched <- match.data(gm.out)
datamatched[datamatched == 99999] <- NA

save(datamatched, file = "datamatched.Rdata")

#__________________________ Table 2, post-matching (genetic)_________________________#

EV.post <- datamatched[2]

covariates.post <- datamatched[, covariate.lbls]

tab2.post <- matrix(NA, nrow = n.covariates, ncol = 4)
rownames(tab2.post) <- covariate.lbls
colnames(tab2.post) <- c("ev", "tv", "diff", "pvalue")

tab2.post[, 1:2] <- cbind(apply(covariates.post[EV.post == 1, ], 2, mean), apply(covariates.post[EV.post == 0,], 2, mean))
tab2.post[, 3] <- tab2.post[, 1] - tab2.post[, 2]
for (i in c(1, 2, 6 , 7)){
  tab2.post[i, 4]<-ks.boot(covariates.post[,i][EV.post==1],covariates.post[,i][EV.post==0], nboots = 500)$ks.boot.pvalue
}
for (i in c(3, 4, 5)){
  tab2.post[i, 4] <- prop.test(table(covariates.post[, i], EV.post$EV), n = apply(table(covariates.post[, i], EV.post$EV),2 , sum))$p.value
}

tab2 <- cbind(tab2.pre, tab2.post)
tab2[3:5, c(1:3, 5:7)] <- tab2[3:5, c(1:3, 5:7)] * 100

### Table 2 ###

print(tab2, digits = 2)


#__________________________ Table 3, pre-matching__________________________#

datamatch[datamatch == 99999] <- NA

outcomes.pre <- datamatch[10:18]

tab3.pre <- matrix(NA,nrow = n.outcomes,ncol = 5)
rownames(tab3.pre) <- outcomes.lbls
colnames(tab3.pre) <- c("N", "prop.ev", "prop.tv", "diff", "pvalue")

for (i in 1:n.outcomes) {
  tab3.pre[i, 1] <- length(na.omit(outcomes.pre[, i]))
  tab3.pre[i, 2:3] <- rev(prop.table(table(outcomes.pre[,i],datamatch$EV),2)[2,])*100
  tab3.pre[i, 4] <- tab3.pre[i, 2] - tab3.pre[i, 3]	
  tab3.pre[i, 5] <- prop.test(table(outcomes.pre[, i], datamatch$EV)[2, ], n = apply(table(outcomes.pre[, i], datamatch$EV), 2, sum))$p.value
}

datamatch[, 10:18][is.na(datamatch[, 10:18]) == "TRUE"] <- 99999


#__________________________ Table 3, post-matching (genetic) _________________________#

outcomes.post <- datamatched[10:18]

tab3.post <- matrix(NA, nrow = n.outcomes, ncol = 5)
rownames(tab3.post) <- outcomes.lbls
colnames(tab3.post) <- c("N", "prop.ev", "prop.tv", "diff", "pvalue")

for (i in 1:n.outcomes) {
  tab3.post[i, 1] <- length(na.omit(outcomes.post[, i]))
  tab3.post[i, 2:3] <- rev(prop.table(table(outcomes.post[, i], datamatched$EV), 2)[2, ]) * 100
  tab3.post[i, 4] <- tab3.post[i, 2] - tab3.post[i, 3]	
  tab3.post[i, 5] <- prop.test(table(outcomes.post[, i], datamatched$EV)[2, ], n = apply(table(outcomes.post[, i], datamatched$EV), 2, sum))$p.value
}

tab3 <- cbind(tab3.pre, tab3.post)

tab3 <- tab3[rev(order(tab3[, 9])), ]

### Table 3 ###

print(tab3, digits = 4)


#______Looking at the balance plot of the propensity score matching________#

library(ggplot2)
library(cobalt)

age_group_propens = bal.plot(m.out, "age.group", which="both")
educ_propens <- bal.plot(m.out, "educ", which="both")
white_collar_propens <- bal.plot(m.out, "white.collar", which="both")
not_full_time_propens <- bal.plot(m.out, "not.full.time", which="both")
male_propens <- bal.plot(m.out, "male", which="both")
tech_propens <- bal.plot(m.out, "tech", which="both")
pol_info_propens <- bal.plot(m.out, "pol.info", which="both")
grid.arrange(age_group_propens, educ_propens, white_collar_propens, not_full_time_propens, male_propens, tech_propens, pol_info_propens, nrow = 4)

#______Looking at the balance plot after genetic matching________#

age_group_propens1 = bal.plot(gm.out, "age.group", which="both")
educ_propens1 <- bal.plot(gm.out, "educ", which="both")
white_collar_propens1 <- bal.plot(gm.out, "white.collar", which="both")
not_full_time_propens1 <- bal.plot(gm.out, "not.full.time", which="both")
male_propens1 <- bal.plot(gm.out, "male", which="both")
tech_propens1 <- bal.plot(gm.out, "tech", which="both")
pol_info_propens1 <- bal.plot(gm.out, "pol.info", which="both")
grid.arrange(age_group_propens1, educ_propens1, white_collar_propens1, not_full_time_propens1, male_propens1, tech_propens1, pol_info_propens1, nrow = 4)

