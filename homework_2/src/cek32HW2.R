###
### Claire Kelling
### STAT 540
### Homework #2
###
### Created 2/21/2018 for the second assigment, due 3/15/2018
### The exercises mainly focus on algorithms for approximating using Bayesian Inference.
### We also use data augmentation and auxiliary variables in this exercise.
### 

library(MASS)
library(mvtnorm)
library(tictoc)
library(optimr)
library(ggplot2)
library(xtable)
library(tictoc)
library(RSpectra)
library(data.table)

###
### Laplace Approximation
###
#clear workspace
rm(list=ls())

#load data, make sure to change to READ FROM WEBSITE *******
#http://personal.psu.edu/muh10/540/Rcode/hw1.dat
prob_5_dat <- fread("http://personal.psu.edu/muh10/540/Rcode/hw1.dat")
y <- c((t(as.matrix(prob_5_dat))))
n <- length(y)

rm(prob_5_dat)

#adopting code from posted solution to homework 1, assuming it is more efficient than mine
#setting variables for laplace
sum.y = sum(y)
log.sum.y = sum(log(y))
alpha.num = function(ab){
  a = ab[1]
  b = ab[2]
  return( -(log(a)+n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y -
              (1/6)*a^2 - (1/6)*b^2 ) )
}
beta.num = function(ab){
  a = ab[1]
  b = ab[2]
  return( -(log(b)+n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y -
              (1/6)*a^2 - (1/6)*b^2 ) )
}
denom = function(ab){
  a = ab[1]
  b = ab[2]
  return( -(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y -
              (1/6)*a^2 - (1/6)*b^2 ) )
}

n.samp <- 10000
a.hat <- rep(NA, n.samp)
b.hat <- rep(NA, n.samp)

tic()
a.num.opt = optim(c(0.01,0.01), alpha.num, hessian=T, method="L-BFGS-B", lower=c(0,0))
b.num.opt = optim(c(0.01,0.01), beta.num, hessian=T, method="L-BFGS-B", lower=c(0,0))
den.opt = optim(c(0.01,0.01), denom, hessian=T, method="L-BFGS-B", lower=c(0,0))

a.num = exp(-alpha.num(a.num.opt$par))*det(a.num.opt$hessian)^{-1/2}
b.num = exp(-beta.num(b.num.opt$par))*det(b.num.opt$hessian)^{-1/2}
den = exp(-denom(den.opt$par))*det(den.opt$hessian)^{-1/2}

a.hat = a.num/den
b.hat = b.num/den
#timing and estimates
lap_time <- toc()
lap_time <- lap_time$toc-lap_time$tic
lap_mean_a <- mean(a.hat)
lap_err_a <- sd(a.hat)/sqrt(n.samp)
lap_mean_b <- mean(b.hat)
lap_err_b <- sd(b.hat)/sqrt(n.samp)


####
#### Problem 1
####

## (Return to the previous homework problem.) Let X1,..,Xn iid from Gamma(a,b) distribution. 
## Consider Bayesian inference for a,b with prior distributions N(0; 3) for both. 
## The data for this problem are available here: http://personal.psu.edu/muh10/540/Rcode/hw1.dat

#(a) Use importance sampling to approximate the expectations. Provide pseudocode, including relevant
# distributions. Is the variability of your Monte Carlo approximation guaranteed to be finite? 
# Explain your answer. Provide Monte Carlo standard errors for your estimates.

n <- 1000
n.rep <- 10
a_init <- 3
b_init <- 0.5
IS.estimates <- rep(NA, n.rep)

#http://dept.stat.lsa.umich.edu/~jasoneg/Stat406/lab7.pdf
for (i in 1:n.rep) {
    # sample size
    n <- length(y)
    # posterior derived in writeup
    log.posterior <- function(alpha) ((n*alpha)*log(beta) - n*log(gamma(alpha)) + sum(log(y^(alpha-1))) -beta*sum(y)+ alpha^2+beta^2)
    # parameters for the trial distribution
    alpha = a_init
    beta = b_init
   
     # log trial density, g
     log.g <- function(t) dgamma(t,alpha,beta,log=TRUE)
     # log importance function
     log.w <- function(t) log.posterior(t) - log.g(t)
     
     # generate from the trial distribution
     res <- 1000
     U <- rgamma(res, alpha, beta)
     # calculate the list of log.w values
     LP <-  log.w(U)
     
     # importance sampling estimate
     IS.estimates[i] <- mean( exp(LP)*U )/mean(exp(LP))
}


tic()
for (i in 1:n.rep) {
  # sample size
  n <- length(y)
  # posterior derived in writeup
  log.posterior <- log(beta^((n*alpha)/(gamma(alpha)^n) + sum(log(y^(alpha-1))) -beta*sum(y)+ alpha^2+beta^2
  # parameters for the trial distribution
  a = a_init; b = b_init
  
  # log trial density, g
  log.g <- function(t) dgamma(t,a,b,log=TRUE)
  # log importance function
  log.w <- function(t) log.posterior(t) - log.g(t)

  # generate from the trial distribution
  res <- 1000
  U <- rgamma(res, a, b)
  # calculate the list of log.w values
  LP <-  log.w(U)

  # importance sampling estimate
  I <- mean( exp(LP)*U )/mean(exp(LP))
  
  # calculate ESS
  ESS <- mean( ( (exp(LP)/mean(exp(LP))) - 1)^2 )
  # calculate s.e. of the IS estimate
  Z = exp(LP)/mean(exp(LP))
  sig.sq <- (1/res)*sum( (Z-I)^2 )
  se <- sqrt( sig.sq/res )
}
imp_time <- toc()
imp_time <- imp_time$toc-imp_time$tic

imp_mean_a <- mean(IS.estimates)
imp_mean_b <- mean(IS.estimates)
imp_err_a <- mcse(IS.estimates, method = "tukey")$se
imp_err_b <- mcse(IS.estimates, method = "tukey")$se

#(b) Construct an all-at-once Metropolis-Hastings (AMH) algorithm to approximate the posterior 
# expectation of a; b. Provide pseudocode for your algorithm. This should include any distributions 
# you had to derive. How did you determine the Monte Carlo approximations were good enough? 
# That is, discuss stopping criteria (how you determined chain length), the number of starting 
# values you tried, how you obtained initial values etc.

proposal=function(x1,x2,v1,v2){
  return(rmvnorm(1, mean = c(x1, x2), sigma=c(v1,v2)*diag(2)))
}

target=function(alpha,beta){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(log(beta^((n*alpha)/(gamma(alpha)^n) + sum(log(y^(alpha-1))) -beta*sum(y)+ alpha^2+beta^2)
  }
}


#Setting intial values and initializing variables
n=1500
mcmc.samp.a=rep(NA,n)
mcmc.samp.b=rep(NA,n)
lambda=3
beta=3
mcmc.samp.a[1]=alpha
mcmc.samp.b[1]=beta

#tuning parameter
set.seed(148752)
v1=0.2
v2=7


for(i in 1:n){
  if(i %% 100 == 0) cat("Starting iteration", i, "\n")
  #proposal=function(x1,x2,v)
  #target=function(lambda,beta)
  x.star=proposal(alpha,beta,v1,v2)
  x1.star=x.star[1]
  x2.star=x.star[2]
  
  log.num=target(x1.star,x2.star)
  log.denom=target(alpha,beta)
  r=exp(log.num-log.denom)
  
  if(runif(1)<r){
    alpha=x1.star
    beta=x2.star
  }
 vmh.mcmc.samp.a[i]=alpha
  mcmc.samp.b[i]=beta
}
amh_time <- toc()
amh_time <- amh_time$toc-amh_time$tic

## Diagnostics
mean(!duplicated(mcmc.samp.a)) # Acceptance rate
mean(!duplicated(mcmc.samp.b))

plot(mcmc.samp.a, type = "l")
plot(mcmc.samp.b, type = "l")

burnin <- 500
mcmc.samp.a <-vmh.mcmc.samp.a[-(1:burnin)]
mcmc.samp.b <- mcmc.samp.b[-(1:burnin)]

plot(mcmc.samp.a, type = "l")
plot(mcmc.samp.b, type = "l")

acf(mcmc.samp.a)
acf(mcmc.samp.b)

library(coda)
effectiveSize(mcmc.samp.a)
effectiveSize(mcmc.samp.b)

amh_mean_a <- mean(mcmc.samp.a)
amh_mean_b <- mean(mcmc.samp.b)
amh_err_a <- mcse(mcmc.samp.a, method = "tukey")$se
amh_err_b <- mcse(mcmc.samp.b, method = "tukey")$se



#(c) Construct a variable-at-a-time Metropolis-Hastings (VMH) algorithm. You need to provide the 
# same level of detail here as you did for the previous algorithm.

# Problem 2
# Part b
proposal1=function(x1,v1){
  return(rnorm(1, mean = x1, sd = sqrt(v1)))
}

proposal2=function(x2,v2){
  return(rnorm(1, mean = x2, sd = sqrt(v2)))
}

target1=function(alpha,beta){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(log(beta^((n*alpha)/(gamma(alpha)^n) + sum(log(y^(alpha-1))) -beta*sum(y)+ alpha^2)
  }
}


target2=function(alpha,beta){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(log(beta^((n*alpha)/(gamma(alpha)^n) + sum(log(y^(alpha-1))) -beta*sum(y)+ beta^2)
  }
}


#Setting intial values and initializing variables
n=1500
vmh.mcmc.samp.a=rep(NA,n)
vmh.mcmc.samp.b=rep(NA,n)
lambda=3
beta=3
vmh.mcmc.samp.a[1]=alpha
vmh.mcmc.samp.b[1]=beta

#tuning parameter
set.seed(89)
#90
v1=3
v2=3
#32- 174 ESS, 42%

tic()
for(i in 1:n){
  if(i %% 100 == 0) cat("Starting iteration", i, "\n")
  #proposal=function(x1,v)
  #target=function(lambda,beta)
  beta.star=proposal1(beta,v)
  log.num=target(alpha,beta.star)
  log.denom=target(alpha,beta)
  r=exp(log.num-log.denom)
  if(runif(1)<r){
    beta=beta.star
  }
  
  alpha.star=proposal1(alpha,v)
  log.num=target(alpha.star,beta)
  log.denom=target(alpha,beta)
  r=exp(log.num-log.denom)
  if(runif(1)<r){
    alpha=alpha.star
  }
  
 vmh.mcmc.samp.a[i]=alpha
 vmh.mcmc.samp.b[i]=beta
}

vmh_time <- toc()
vhm_time <- vmh_time$toc-vmh_time$tic

## Diagnostics
# Acceptance rate
mean(!duplicated(vmh.mcmc.samp.a))
mean(!duplicated(vmh.mcmc.samp.b))

plot(vmh.mcmc.samp.a, type = "l")
plot(vmh.mcmc.samp.b, type = "l")

burnin <- 500
vmh.mcmc.samp.a <-vmh.mcmc.samp.a[-(1:burnin)]
vmh.mcmc.samp.b <- mcmc.samp.b[-(1:burnin)]

plot(vmh.mcmc.samp.a, type = "l")
plot(vmh.mcmc.samp.b, type = "l")

acf(vmh.mcmc.samp.a)
acf(vmh.mcmc.samp.b)

library(coda)
effectiveSize(vmh.mcmc.samp.a)
effectiveSize(vmh.mcmc.samp.b)

vmh_mean_a <- mean(vmh.mcmc.samp.a)
vmh_mean_b <- mean(vmh.mcmc.samp.b)
vmh_err_a <- mcse(vmh.mcmc.samp.a, method = "tukey")$se
vmh_err_b <- mcse(vmh.mcmc.samp.b, method = "tukey")$se


#(d) Provide a table with all approximations along with any error approximations, as well as the 
# computational time taken by the algorithms. (To make this easy to read, the rows of your table 
# should correspond to algorithms.)
lap_col <- c(lap_mean_a, lap_mean_b, lap_err_a, lap_err_b, lap_time)
imp_col <- c(imp_mean_a, imp_mean_b, imp_err_a, imp_err_b, imp_time)
amh_col <- c(amh_mean_a, amh_mean_b, amh_err_a, amh_err_b, amh_time)
vmh_col <- c(vmh_mean_a, vmh_mean_b, vmh_err_a, vmh_err_b, vmh_time)

full_df <- cbind(lap_col, imp_col, amh_col, vmh_col)
rownames(full_df) <- c("a_est", "b_est", "a_err", "b_err", "approx_time")
colnames(full_df) <- c("Laplace", "Imp Samp", "AMH", "VMH")


#(e) Compare the efficiency of the four algorithms {importance sampling, AMH, VMH, and the
# Laplace approximation from the previous homework} using the methodology discussed in class, 
# e.g. effective sample size, effective samples per second etc. The Monte Carlo algorithms are 
# easier to compare than the comparison between them and the Laplace approximation; you may 
# have to think carefully about this. Be sure to clearly explain why your approaches for comparing
# the algorithms are reasonable.

#Laplace

#IMP
effectiveSize(IS.estimates.a)
effectiveSize(IS.estimates.b)


#AMH
effectiveSize(mcmc.samp.a)
effectiveSize(mcmc.samp.b)

#VMH
effectiveSize(vmh.mcmc.samp.a)
effectiveSize(vmh.mcmc.samp.b)

#(f) Which algorithm would you recommend for this problem? How would you order the algorithms
# in terms of ease of implemention?



###
### Problem 2
###

# Lightbulbs: Assume that lightbulb lifetimes for a lightbulb made by a particular company are 
# independent and exponentially distributed with expectation theta. Suppose in an experiment, 
# m bulbs are switched on at the same time, but are only completely observed up to time tau .
# Let the lifetimes of these bulbs be A1, A2, ..., Am. However, since the bulbs are only 
# observed till time tau , not all these lifetimes will be observed. Now suppose that at time 
# tau , the experimenter observes the number of lightbulbs, W, still working at time tau, and 
# the lifetimes of all lightbulbs that stopped working by tau . For convenience, denote these 
# bulbs that stopped working by time tau as A1, A2,...,Am????W. Hence, the missing informa  to
# n consists of the lifetimes of the lightbulbs still working at time tau, Am????W+1; : : : ;Am.
# For a particular experiment, let tau be 184 days and m = 100. The observations for this 
# experiment are as follows. The data on the lightbulb lifetimes for the bulbs that stopped 
# working by tau are on the website. 
# 
# (Assume that the remaining bulbs were still working at time tau.) Let the prior for theta be
# Gamma(1,100) (with parameterization so it has mean 100 and variance 100^2).

#load data
bulb_dat <- fread("http://personal.psu.edu/muh10/540/hwdir/bulbs.dat")
bulb_dat <- c((t(as.matrix(bulb_dat))))

#(a) Using auxiliary variables for the missing lightbulb data, construct an MCMC algorithm 
# to approximate the posterior distribution of theta. Provide the same level of detail about 
# your MCMC algorithm as you did in the previous problem.

#prior for theta is Gamma(1,100)
#lifetime of lightbulbs are iid exp(theta)
#let tau be 184 days and m = 100
times <- c((t(as.matrix(bulb_dat))))

for( i in 1:n.rep){
  
  
}


#b) Construct a different MCMC algorithm, this time without using auxiiliary variables/data 
# augmentation. Again, provide details.

# This time, I will approximate the posterior directly, without data augmentation.
 
#prior for theta is Gamma(1,100)
#lifetime of lightbulbs are iid exp(theta)
#let tau be 184 days and m = 100
#posterior distribution of theta

n.rep <- 1000
n <- length(bulb_dat)

mcmc_wo_aux <- rgamma(n.rep, shape = n+1, scale = (sum(bulb_dat)+ 1/100)^(-1))
hist(mcmc_wo_aux)
mean(mcmc_wo_aux)


#(c) Overlay posterior density plot approximations for the two algorithms. Provide a table 
# that shows the posterior mean approximations for tau along with MCMC standard errors.
plot(density(mcmc_w_aux), "Density with auxiliarly variables")
lines(density(mcmc_wo_aux))

mean_aux <- mean(mcmc_w_aux])
mean_wo_aux <- mean(mcmc_wo_aux)
se_aux <- sd(mcmc_w_aux)/sqrt(n.samp)
se_wo_aux <- sd(mcmc_wo_aux)/sqrt(n.samp)

#(d) For the auxiliary variable method, plot the approximate posterior pdf for one of the 
# "missing" lightbulbs, then overlay the approximate posterior pdfs for a lightbulb made 
# by the company. You should notice that they are different. Report the posterior mean 
# estimates for each of them. 


#(e) Looking at just the ability to approximate the posterior per iteration of the algorithm, 
# which of the two MCMC algorithms is more efficient? Now accounting for computing costs, 
# which of the two MCMC algorithms is more efficient? Which algorithm would you recommend?



#(f) This course is focused on computing but it is worth noting some basics about inference. 
# Compare your results above with what would happen to inference if you ignored the missing 
# data by overlaying the density plots.
