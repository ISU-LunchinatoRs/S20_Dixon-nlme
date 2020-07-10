# nl me example - data from Fernando's nlraa package

# live fuel moisture content
#   only for one species: S. bracteolactus
#   3 plots, data measured over time, 7 times

library(dplyr)    # some data manipulation
library(lattice)  # data plots by plot

library(nlraa)    # source of data
library(nlme)     # nlme()

library(lme4)     # nlmer()
library(metafor)  # meta analysis
library(rstanarm) # Bayesian nlmer

data(lfmc)
sapply(lfmc, class)

with(lfmc, table(site, plot))

# consider only S. bracteolactus
sb <- lfmc %>% filter(leaf.type=="S. bracteolactus")
with(sb, table(plot, time))

xyplot(lfmc ~ time | plot, data=sb, layout=c(3,1))

# using nlme functions

# fit fixed effects models to each plot

# first define the y, x, and grouping variables
#   group must be a factor variable (already set up)
#   can specify nested groups (B within A) by A/B
#     (more later)
#   result depends on whether X is continuous or factor
sbGrp <- groupedData(lfmc ~ time | plot, data=sb)

# and there is a nice lattice graph plot of data in each group
plot(sbGrp)

# get all individual fits with nlsList
# SSxxxx functions are "self-start" nonlinear models, include 
#   usually reasonable estimates of starting values
# Fernando's vignette for nlraa lists the SS functions in nlme and nlraa
#   can write your own
sb.plot <- nlsList(lfmc ~ SSdlf(time, upper, lower, mid, scale), data=sbGrp)
sb.plot
summary(sb.plot)
intervals(sb.plot)

par(mfrow=c(1,3), mar=c(3,3,0,0)+0.3, mgp=c(2,0.8,0))
sb.pred1 <- predict(sb.plot)
for (i in unique(sbGrp$plot)) {
  keep <- sbGrp$plot==i
  plot(sbGrp$time[keep], sbGrp$lfmc[keep], pch=19, col=4)
  lines(sbGrp$time[keep], sb.pred1[keep])
}


# predicting every day to get smoother plot
# requires a bit more setup
days <- expand.grid(time=0:80, plot=1:3)
# include plot so can extract for each plot
days$pred <- predict(sb.plot, newdata=days)
# and store prediction in that data frame
# one quirk - plot variable ignored; group info used instead
#   so days needs to have all times for plot 1, then for plot 2
#   order of variables in the expand.grid() matters

# if have
#  days <- expand.grid(plot=1:3, time=0:80)
# you clearly get the wrong fitted curve

for (i in unique(sbGrp$plot)) {
  keep <- sbGrp$plot==i
  plot(sbGrp$time[keep], sbGrp$lfmc[keep], pch=19, col=4)
  keeppred <- days$plot==i
  lines(days$time[keeppred], days$pred[keeppred])
}

# fit nl mixed effect model
# easy if start with the plot-specific fits (the nlsList object)

# get a warning message - don't ignore
sb.nlme <- nlme(sb.plot)

# blindly follow the advice and increase # iterations

sb.nlme <- nlme(sb.plot, 
  control=nlmeControl(MaxIter = 200) )
sb.nlme <- nlme(sb.plot, 
  control=nlmeControl(MaxIter = 200, msMaxIter=200) )

# instead of blindly trying harder, let's look at the output:
sb.nlme

# we see that nlme is trying to fit 4 random effects (one per parameter)
#   with arbitrary correlation matrix 
#   General positive-definite or look at estimated RE structure: has correlations
# trying to estimate 4 variances and 6 correlations from 3 levels (plots)
# this is the default

# simplify the random effects structure - no correlations
# that is a pdDiag() variance covariance matrix
sb.nlme2 <- nlme(sb.plot, random = pdDiag(upper + lower + mid + scale ~ 1) )

# this runs without error
sb.nlme2

# could also use update(sp.nlme, random = pdDiag(upper + lower + mid + scale ~ 1) )
#   to change the random effect specification without specifying everything again

# ?pdClasses tells you the various possibilities
#    default is pdLogChol, which is a better way to parameterize pdSymm

# look at the output
summary(sb.nlme2)

# confidence intervals for the fixed effects
intervals(sb.nlme2, which='fixed')

# look at the plot-specific coefficients: two ways

# estimated fixed effects and predicted random effects
fixef(sb.nlme2)
ranef(sb.nlme2)

# matrix of coefficients, rows = plots
coef(sb.nlme2)

# look at residual vs predicted value plot
# these are standardized (variance = 1) conditional residuals 
#   i.e. given BLUPs of the random effects
# useful to validate assumptions about the error distribution
plot(sb.nlme2)

# plot the results: two possibilities:
#  predict using the fixed effects (same curve for all three plots)
plot(augPred(sb.nlme2, level=0))

# or using the plot-specific coefficients
plot(augPred(sb.nlme2, level=1))

# or combine the two
plot(augPred(sb.nlme2, level=0:1))

# results suggest no variability in lower, mid and scale
#   just in upper

# do not need to include all random effects
#   also illustrating starting with a data frame, not the nlsList
# possible with update because model specified in sb.nlme2
sb.nlme3 <- update(sb.nlme2, random = pdDiag(upper ~ 1), data=sbGrp)

# if start with Grouped data, need to include the model
sb.nlme3 <- nlme(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
  random = pdDiag(upper ~ 1), data=sbGrp)

# but you need to provide more information if try to start with data
#   I've always struggled to correctly specify all the pieces.
# my recommendation: start with a grouped data object

# fit model with no random effects using nls()
sb.nlme4 <- nls(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
   data=sbGrp)

# which random effect model is more appropriate?
c(corr = AIC(sb.nlme), var4 = AIC(sb.nlme2), 
  var1=AIC(sb.nlme3), fixed=AIC(sb.nlme4) )
 c(corr = BIC(sb.nlme), var4 = BIC(sb.nlme2), 
  var1=BIC(sb.nlme3), fixed=BIC(sb.nlme4) ) 
 
# what if there are multiple levels of nesting

sb2 <- read.csv('lfmcID.csv', as.is=T)
sb2$plot <- factor(sb2$plot)
sb2$ID <- factor(sb2$ID)

sb2.Grp <-  groupedData(lfmc ~ time | plot/ID, data=sb2)

# models fit separately to each plot and ID
sb2.plotID <- nlsList(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
  data=sb2.Grp)
sb2.plotID

# fitting a nl me with nested groups: code from Fernando
#   variability only in upper

# Starting with the grouped data doesn't work
sb2.nlme <- nlme(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
  random = list(plot = upper ~ 1, ID= upper~1), 
  groups = ~ plot/ID,
  data=sb2.Grp)

# starting with a single effect model, then updating does
sb1.nlme <- nlme(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
  random = upper ~ 1,
  data=sb2.Grp)

# but you have to turn off the "nested groups" first
sb2.GrpB <-  groupedData(lfmc ~ time | plot, data=sb2)

sb1.nlme <- nlme(lfmc ~ SSdlf(time, upper, lower, mid, scale), 
  random = upper ~ 1,
  data=sb2.GrpB)

sb2.nlme <- update(sb1.nlme, 
  random = list(group = upper ~ 1,
     ID = upper ~ 1),
  groups = ~ group/ID)

sb2.nlme

# and can extend to different RE structure at each level
sb3.nlme <- update(sb1.nlme, 
  random = list(group = pdMat(upper + scale ~ 1),
     ID = pdMat(upper ~ 1) ),
  groups = ~ group/ID)

sb3.nlme

# if the normal approximation (for estimates) is suspect, use a parametric bootstrap

# nlme includes a simulate.lme() function, but that fails with a non-linear fit
sb.sim <- simulate(sb.nlme3, method='ML')

# can also do much of this with nlmer
# differences: 
# 1) NL function must return predictions as a vector 
#    AND have a gradient attribute.  Self-start functions include gradients
# 1) specify random effects in a 3 part formula
#   response ~ NL model ~ 
# 2) no need to group data first
# 3) self-starting part of SSxxxx functions not (currently) used
#   need to specify starting values as a named vector
#   fixed effect estimates ignoring plot is often a good start

sb.nlmer3 <- nlmer(
  lfmc ~ SSdlf(time, upper, lower, mid, scale) ~ upper | plot,
  start=c(upper=286, lower=53, mid=33, scale=-16),
  data=sb)
# complains about not converged, change optimizer
# Default is Nelder-Mead, which is often slow and troublesome
# if needed more interations, add , optCtrl=list(maxfun=20000)
#   to the nlmerControl argument

sb.nlmer3 <- nlmer(
  lfmc ~ SSdlf(time, upper, lower, mid, scale) ~ upper | plot,
  start=c(upper=286, lower=53, mid=33, scale=-16),
  control=nlmerControl(optimizer='bobyqa'),
  data=sb)

# could also use the SSfpl(): four param logistic function from base stats
sb.nlmer3b <- nlmer(
  lfmc ~ SSdlf(time, upper, lower, mid, scale) ~ upper | plot,
  start=c(upper=286, lower=53, mid=33, scale=-16),
  control=nlmerControl(optimizer='bobyqa'),
  data=sb)

# can add additional random effects, either correlated or not
# correlated (bad here, because corr = 1)
sb.nlmer4a <- nlmer(
  lfmc ~ SSdlf(time, upper, lower, mid, scale) ~ (upper + lower | plot),
  start=c(upper=286, lower=53, mid=33, scale=-16),
  control=nlmerControl(optimizer='bobyqa'),
  data=sb)

sb.nlmer4b <- nlmer(
  lfmc ~ SSdlf(time, upper, lower, mid, scale) ~ 
    (upper | plot)  + (lower | plot),
  start=c(upper=286, lower=53, mid=33, scale=-16),
  control=nlmerControl(optimizer='bobyqa'),
  data=sb)

# lme4 includes a bootMer() function for bootstrap 
#   confidence intervals
# I haven't yet figured out how to use it correctly

# comparisons of nlme and nlmer predictions

# for plots in the data set
# predictions - for plots in the data set
#   use blups of the random effects for each plot

sb.nlme.pred <- predict(sb.nlme3)
# default is finest grouping level

sb.nlmer.pred <- predict(sb.nlmer3)

# compare them
par(mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
par(mfrow=c(1,3))
for (i in unique(sb$plot)) {
  bit <- subset(sb, plot==i)
  plot(bit$time, bit$lfmc, pch=19, col=4, 
    xlab='Time', ylab='lfmc')
  lines(bit$time, sb.nlme.pred[sb$plot==i], lty=1, col=3, lwd=2)
  lines(bit$time, sb.nlmer.pred[sb$plot==i], lty=2, col=4, lwd=2)
  legend('topright',bty='n', lty=1:2, col=c(4,3), 
    legend=c('nlme', 'nlmer'), lwd=2)
  }

# predictions for a new plot
#   based on fixed effect model

newdata <- data.frame(time=1:72)

sb.nlme.pred0 <- predict(sb.nlme3, newdata=newdata, level=0)
# level 0 is the population estimates of fixed effects

# Should be possible for predictions from nlmer objects
#   but I haven't figured it out (yet)

# Alternatives to frequentist mixed effect models

# meta analysis:
#   requires that can fit to each subject (e.g. plot)
#   starts with plot-specific estimates and standard errors

# easiest way (that I know) is to extract from the 
#  coefficients part of the summary

temp <- summary(sb.plot)$coef
temp

# first column is the estimate, second is the se
upper <- temp[,1:2, 'upper']
upper

# random effects meta analysis
rma(yi=upper[,1], sei=upper[,2], method='REML')

# fixed effects meta analysis
rma(yi=upper[,1], sei=upper[,2], method='FE')

# differences from nlme:
# 1) MA estimates based on plot-specific model
#     all parameters differ among plots
#   nlme can fit models with some parameters differing
#    others same for all plots (omitted from random = )
# 2) MA uses only moments (estimates, se's)
#   nlme uses full distribution of random effects
#     experience => not very sensitive to non-normal re's

# Bayesian inference using rstanarm
# can ONLY use self starting functions provided in base stats:
#    SSasymp, SSasympOff, SSasympOrig, SSbiexp, SSfol, SSfpl, SSgompertz, 
#    SSlogis, SSmicmen, and SSweibull.
sb.bayes <- stan_nlmer(
  lfmc ~ SSfpl(time, upper, lower, mid, scale) ~ upper | plot,
#  control=nlmerControl(optimizer='bobyqa'),
  chains=3,
  iter=5000,
  data=sb)

# check convergence using Rhat
summary(sb.bayes)[,'Rhat']

# lots of information about estimates and diagnostics
summary(sb.bayes)

# visual exploration of model results 
#  if a large data set, probably want to turn off
# posterior predictive checks (   ,ppd = F)
launch_shinystan(sb.bayes)

# extract plot-specific coefficients
coefficients(sb.bayes)

posterior_interval(sb.bayes)

pairs(sb.bayes, 
  pars=c('upper','lower','mid','scale', 
         'Sigma[plot:upper,upper]'))

# posterior distributions are far from normal
#  Wald inference (i.e. nlme or nlmer is suspect)
#  profile likelihood intervals would be ok
#    but not yet available in nlmer


