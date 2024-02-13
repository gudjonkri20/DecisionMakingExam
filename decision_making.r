seed_id = 1982
set.seed(seed_id)

install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

rawDat <- read.csv("C:\\Users\\gudjo\\Downloads\\HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

#- create covariates in raw data matrix
# nation index
rawDat$nation <- c()
rawDat$nation[rawDat$city=="Melbourne"]=1
rawDat$nation[rawDat$city=="Minsk"]=2
rawDat$nation[rawDat$city=="Chengdu"]=3
rawDat$nation[rawDat$city=="Copenhagen"]=4
rawDat$nation[rawDat$city=="Bonn"]=5
rawDat$nation[rawDat$city=="Athens"]=6
rawDat$nation[rawDat$city=="Seoul"]=7
rawDat$nation[rawDat$city=="Samara"]=8
rawDat$nation[rawDat$city=="Muscat"]=9
rawDat$nation[rawDat$city=="Riyadh"]=10
rawDat$nation[rawDat$city=="Istanbul"]=11
rawDat$nation[rawDat$city=="Nottingham"]=12
rawDat$nation[rawDat$city=="Dnipropetrovs'k"]=13
rawDat$nation[rawDat$city=="Boston"]=14

# create variable for GINI. Data from 
# http://hdr.undp.org/sites/default/files/reports/269/hdr_2009_en_complete.pdf,

rawDat$gini <- c()
rawDat$gini[rawDat$city=="Melbourne"]=12.3
rawDat$gini[rawDat$city=="Minsk"]=1.5
rawDat$gini[rawDat$city=="Chengdu"]=11.6
rawDat$gini[rawDat$city=="Copenhagen"]=26.2
rawDat$gini[rawDat$city=="Bonn"]=9.3
rawDat$gini[rawDat$city=="Athens"]=64.0
rawDat$gini[rawDat$city=="Seoul"]=6.7
rawDat$gini[rawDat$city=="Samara"]=12.1
rawDat$gini[rawDat$city=="Muscat"]=44.7
rawDat$gini[rawDat$city=="Riyadh"]=36.7
rawDat$gini[rawDat$city=="Istanbul"]=17.7
rawDat$gini[rawDat$city=="Nottingham"]=40.5
rawDat$gini[rawDat$city=="Dnipropetrovs'k"]=8.8
rawDat$gini[rawDat$city=="Boston"]=6.3

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)
singlechild <- matrix(NA, nrow = groupSize, ncol = ngroups)


for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    singlechild[s, g] <- redDat$singlechild[redDat$groupid == group_names[g] & redDat$p == "N-experiment"][s]
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

Gini <- array(0, ngroups)
Nation <- array(0, ngroups)
for (g in 1:ngroups) {
  Gini[g] <- mean(redDat$gini[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

c_win <- c_no_punish[,,!is.na(Gini)]
c_keep <- rep()-c_win

Gga_punish <- Gga_punish[,!is.na(Gini)]
Gga_no_punish <- Gga_no_punish[,!is.na(Gini)]

c <- c[,,!is.na(Gini),]
Gga <- Gga[,!is.na(Gini),]
Ggas <- Ggas[,,!is.na(Gini),]
Gc <- Gc[,,!is.na(Gini),]
Ga <- Ga[,,!is.na(Gini),]
Gini <- Gini[!is.na(Gini)]
Nation <- Nation[!is.na(Nation)]

#redefine number of groups after removing those without civic scores
ngroups <- length(Gini)

# aggregate Gini to just 1 number per Nation-index (using the mean here should be unproblematic since all groups within a given nation should have been given the same Gini-coefficient)
Gini <- aggregate(Gini~Nation, FUN=mean)[,2]

nnations <- length(Gini)

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

#################################################################
#------------------ Winnings analysis ---------------------------
#################################################################

#-------------------  Regress Gini on winnings ---------------

# standardise variables

X <- Gini
X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

Y <- (winnings-mean(winnings))/sd(winnings)

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# - run jags code
win.samples <- jags.parallel(data, inits=NULL, params,
                    model.file ="C:\\Users\\gudjo\\OneDrive\\Desktop\\Cognitive Science\\Decision Making\\win_corr.txt",
                    n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  Regress Gini on belief weights and slope of prefs in CC model ---------------

# standardise covariate
X <- Gini
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

# Ga_old <- Ga
# Ga <- Ggas

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                   model.file ="C:\\Users\\gudjo\\OneDrive\\Desktop\\Cognitive Science\\Decision Making\\CC_corr.txt",
                   n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)
end_time = Sys.time()
end_time - start_time

#################################################################
#------------------ Plotting ---------------------------
#################################################################

# ------ Create empirical and parameter arrays for plots -------
# empirical group winnings data - means and standard deviations
empirical.win <- array(0,c(3,length(Gini)))
for (i in 1:length(Gini)) {
  empirical.win[1,i] <- mean(winnings[Nation==i]) - sd(winnings[Nation==i]) 
  empirical.win[2,i] <- mean(winnings[Nation==i]) 
  empirical.win[3,i] <- mean(winnings[Nation==i]) + sd(winnings[Nation==i])
}

empirical.initial <- array(0,c(3,length(Gini)))
for (i in 1:length(Gini)) {
  empirical.initial[1,i] <- mean(colMeans(c[,1,,1])[Nation==i]) - sd(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[2,i] <- mean(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[3,i] <- mean(colMeans(c[,1,,1])[Nation==i]) + sd(colMeans(c[,1,,1])[Nation==i])
}

#-------------- create high res png for saving ---------------------------------

png(filename="PGG_offshore_tax_figure_1.png", height = 800, width = 1000)
layout(matrix(1:2, nrow = 1))

# Graph A: Data - Contribution
par(mar=c(5,4,4,2) + 0.1)
plot(c(1,65), c(0,1300), type = "n", main = "A: Group Contribution", 
     xlab = "% of National GDP held in tax havens by households", ylab = "Group Contribution", axes=FALSE)
axis(1)
axis(2)
# Replace with your actual plotting code for A
# Assuming 'Gini' and 'empirical.win' are your data
for (i in 1:length(Gini)) {
  lines(c(Gini[i], Gini[i]), c(empirical.win[1,i], empirical.win[3,i]))
  points(Gini[i], empirical.win[2,i], pch=19)
}

# Graph B: Posterior of effect - Gini and group winnings with adjusted 95% interval
par(mar=c(5,4,4,2) + 0.1)
plot(density(win.samples$BUGSoutput$sims.list$betaX), frame=FALSE, lwd=2, ylim=c(0,7),
     main="B: Mean Posterior Contribution", xlab = "Standardised % of National GDP held in tax havens by households", ylab = "")
abline(v=0, lty=2) # Zero line for reference
# Adjusted 95% Bayesian Credible Interval for B
y_pos <- -0.05 # Adjust this value as needed to place the line within visible range but at the bottom
COL <- "red" # Directly use the color
lines(c(win.samples$BUGSoutput$summary[2, "2.5%"], win.samples$BUGSoutput$summary[2, "97.5%"]),
c(y_pos, y_pos), col=COL, lwd=2)
points(win.samples$BUGSoutput$summary[2, "50%"], y_pos, pch=19, col=COL) # Median point
  
dev.off()


png(filename="PGG_offshore_tax_figure_2.png", height = 800, width = 1000)
layout(matrix(1:2, nrow = 1))

# Graph C: Data - Initial Contribution
par(mar=c(5,4,4,2) + 0.1)
plot(c(1,65), c(0,30), type = "n", main = "C: Initial Contribution", 
     xlab = "% of National GDP held in tax havens by households", ylab = "Initial Contribution", axes=FALSE)
axis(1)
axis(2)
# Replace with your actual plotting code for C
# Assuming 'Gini' and 'empirical.initial' are your data
for (i in 1:length(Gini)) {
  lines(c(Gini[i], Gini[i]), c(empirical.initial[1,i], empirical.initial[3,i]))
  points(Gini[i], empirical.initial[2,i], pch=19)
}

# Graph D: Initial Belief with adjusted 95% interval
par(mar=c(5,4,4,2) + 0.1)
plot(density(CC.samples$BUGSoutput$sims.list$betaX_alpha), frame=FALSE, lwd=2, ylim=c(0,15),
     main="D: Initial Belief", xlab = "Standardised % of National GDP held in tax havens by households", ylab = "")
abline(v=0, lty=2) # Zero line for reference
y_pos_d <- -0.5 # Adjust this value as needed
# Adjusted 95% Bayesian Credible Interval for D
lines(c(CC.samples$BUGSoutput$summary[4, "2.5%"], CC.samples$BUGSoutput$summary[4, "97.5%"]),
c(y_pos_d, y_pos_d), col="red", lwd=2)
points(CC.samples$BUGSoutput$summary[4, "50%"], y_pos_d, pch=19, col="red")


dev.off()



png(filename="PGG_offshore_tax_figure_3.png", height = 800, width = 1000)
layout(matrix(1:2, nrow = 1))

# Graph E: Belief Learning Weight with adjusted 95% interval
par(mar=c(5,4,4,2) + 0.1)
plot(density(CC.samples$BUGSoutput$sims.list$betaX_omega), frame=FALSE, lwd=2, ylim=c(0,10),
     main="E: Belief Learning Weight", xlab = "Standardised % of National GDP held in tax havens by households", ylab = "")
abline(v=0, lty=2) # Zero line for reference
y_pos_e <- -0.05
# Drawing the 95% Bayesian Credible Interval at the bottom
lines(c(CC.samples$BUGSoutput$summary[5, "2.5%"], CC.samples$BUGSoutput$summary[5, "97.5%"]),
      c(y_pos_e, y_pos_e), col="red", lwd=2)
points(CC.samples$BUGSoutput$summary[5, "50%"], y_pos_e, pch=19, col="red")

# Graph F: Conditional Preferences with adjusted 95% interval
par(mar=c(5,4,4,2) + 0.1)
plot(density(CC.samples$BUGSoutput$sims.list$betaX_rho), frame=FALSE, lwd=2, ylim=c(0,10),
     main="F: Conditional Preferences", xlab = "Standardised % of National GDP held in tax havens by households", ylab = "")
abline(v=0, lty=2) # Zero line for reference
y_pos_f <- -0.05
# Drawing the 95% Bayesian Credible Interval at the bottom
lines(c(CC.samples$BUGSoutput$summary[6, "2.5%"], CC.samples$BUGSoutput$summary[6, "97.5%"]),
      c(y_pos_f, y_pos_f), col="red", lwd=2)
points(CC.samples$BUGSoutput$summary[6, "50%"], y_pos_f, pch=19, col="red")

dev.off()



