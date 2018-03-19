###
### Claire Kelling
### STAT 540
### Homework #2
###
### Created 2/21/2018 for the second assigment, due 3/16/2018
### The exercises mainly focus on algorithms for approximating posteriors using importance sampling,
### all-at-once Metropolis Hastings and variable-at-a-time Metropolis Hastings.
### We also use data augmentation and auxiliary variables in this exercise.
### 

library(MASS)
library(coda)
library(mvtnorm)
library(tictoc)
library(optimr)
library(ggplot2)
library(xtable)
library(tolerance)
library(tictoc)
library(RSpectra)
library(mcmcse)
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

n.samp <- 100000

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


####
#### Problem 1
####

## (Return to the previous homework problem.) Let X1,..,Xn iid from Gamma(a,b) distribution. 
## Consider Bayesian inference for a,b with prior distributions N(0; 3) for both. 
## The data for this problem are available here: http://personal.psu.edu/muh10/540/Rcode/hw1.dat

#(a) Use importance sampling to approximate the expectations. Provide pseudocode, including relevant
# distributions. Is the variability of your Monte Carlo approximation guaranteed to be finite? 
# Explain your answer. Provide Monte Carlo standard errors for your estimates.

IS.estimates.a <- rep(NA, n.samp)
IS.estimates.b <- rep(NA, n.samp)
IS.estimates.a[1] <- a.hat
IS.estimates.b[1] <- b.hat

tic()
for (i in 2:n.samp) {
  if(i %% 1000 == 0) cat("Starting iteration", i, "\n")
  ####
  #### First for alpha
  ####
    # posterior derived in writeup
    log.posterior <- function(a,b){
      return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - (1/6)*a^2 )
      #return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y - (1/6)*a^2 - (1/6)*b^2)
    }

     # log proposal density, g
     log.g <- function(t) dexp(t,rate=1/t,log=TRUE)
     
     # log importance function
     log.w <- function(t) log.posterior(t, IS.estimates.b[i-1]) - log.g(t)
     
     # generate from the proposal distribution
     res <- 1000
     U <- rexp(res, rate = 1/IS.estimates.a[i-1])
     # calculate the list of log.w values
     LP <-  log.w(U)
     
     # importance sampling estimate
     IS.estimates.a[i] <- mean( exp(LP)*U )/mean(exp(LP))
     
     
     ####
     #### Now for beta
     ####
     # posterior derived in writeup
     log.posterior <- function(a,b){
       return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y - (1/6)*a^2 - (1/6)*b^2)
     }

     # log proposal density, g
     log.g <- function(t) dexp(t,rate=1/IS.estimates.b[i-1],log=TRUE)
     
     # log importance function
     log.w <- function(t) log.posterior(IS.estimates.a[i-1],t) - log.g(t)
     
     # generate from the proposal distribution
     res <- 1000
     U <- rexp(res, 1/IS.estimates.b[i-1])
     # calculate the list of log.w values
     LP <-  log.w(U)
     
     # importance sampling estimate
     IS.estimates.b[i] <- mean( exp(LP)*U )/mean(exp(LP))
}

imp_time <- toc()
imp_time <- imp_time$toc-imp_time$tic

#Diagnostics
# mean(!duplicated(IS.estimates.a)) # Acceptance rate, not relevant
# mean(!duplicated(IS.estimates.b))

plot(IS.estimates.a, type = "l", main = "Trace plot for alpha")
plot(IS.estimates.b, type = "l", main = "Trace plot for beta")

burnin <- 500
IS.estimates.a <- IS.estimates.a[-(1:burnin)]
IS.estimates.b <- IS.estimates.b[-(1:burnin)]

plot(IS.estimates.a, type = "l",  main = "Trace plot for alpha, without burnin")
plot(IS.estimates.b, type = "l",  main = "Trace plot for beta, without burnin")

acf(IS.estimates.a, main = "ACF for alpha, without burnin")
acf(IS.estimates.b, main = "ACF for beta, without burnin")

imp_ess_a <- effectiveSize(IS.estimates.a)
imp_ess_b <- effectiveSize(IS.estimates.b)

imp_ess_a_t <- effectiveSize(IS.estimates.a)/imp_time
imp_ess_b_t <- effectiveSize(IS.estimates.b)/imp_time

imp_mean_a <- mean(IS.estimates.a)
imp_mean_b <- mean(IS.estimates.b)
imp_err_a <- mcse(IS.estimates.a, method = "tukey")$se
imp_err_b <- mcse(IS.estimates.b, method = "tukey")$se

df <- rbind(c(imp_mean_a, imp_mean_b), c(imp_err_a, imp_err_b))
xtable(df, digits = c(1,7,7))

#(b) Construct an all-at-once Metropolis-Hastings (AMH) algorithm to approximate the posterior 
# expectation of a; b. Provide pseudocode for your algorithm. This should include any distributions 
# you had to derive. How did you determine the Monte Carlo approximations were good enough? 
# That is, discuss stopping criteria (how you determined chain length), the number of starting 
# values you tried, how you obtained initial values etc.

proposal=function(a,b){
  #indep exponentials- only using same support as alpha,beta
  a_prop <- rexp(1,rate=1/a)
  b_prop <- rexp(1,rate=1/b)
  prop <- c(a_prop, b_prop)
  return(prop)
}

target=function(a,b){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y - (1/6)*a^2 - (1/6)*b^2 )
  }
}


#Setting intial values and initializing variables
mcmc.samp.a=rep(NA,n.samp)
mcmc.samp.b=rep(NA,n.samp)
alpha=a.hat
beta=b.hat
mcmc.samp.a[1]=alpha
mcmc.samp.b[1]=beta

tic()
for(i in 1:n.samp){
  if(i %% 100 == 0) cat("Starting iteration", i, "\n")

  x.star=proposal(alpha,beta)
  x1.star=x.star[1]
  x2.star=x.star[2]
  
  log.num=target(x1.star,x2.star)
  log.denom=target(alpha,beta)
  r=exp(log.num-log.denom)
  
  if(runif(1) < r){
    alpha=x1.star
    beta=x2.star
  }
  mcmc.samp.a[i]=alpha
  mcmc.samp.b[i]=beta
}
amh_time <- toc()
amh_time <- amh_time$toc-amh_time$tic

## Diagnostics
mean(!duplicated(mcmc.samp.a)) # Acceptance rate
mean(!duplicated(mcmc.samp.b))

plot(mcmc.samp.a, type = "l", main ="Trace plots of alpha with burnin")
plot(mcmc.samp.b, type = "l", main ="Trace plots of beta with burnin")

burnin <- 5000
mcmc.samp.a <- mcmc.samp.a[-(1:burnin)]
mcmc.samp.b <- mcmc.samp.b[-(1:burnin)]

plot(mcmc.samp.a, type = "l", main = "Trace plots of alpha without burnin")
plot(mcmc.samp.b, type = "l", main = "Trace plots of beta without burnin")

acf(mcmc.samp.a, main = "ACF of alpha without burnin")
acf(mcmc.samp.b, main = "ACF of beta without burnin")

amh_ess_a <- effectiveSize(mcmc.samp.a)
amh_ess_b <- effectiveSize(mcmc.samp.b)

amh_ess_a_t <- effectiveSize(mcmc.samp.a)/amh_time
amh_ess_b_t <- effectiveSize(mcmc.samp.b)/amh_time

amh_mean_a <- mean(mcmc.samp.a)
amh_mean_b <- mean(mcmc.samp.b)
amh_err_a <- mcse(mcmc.samp.a, method = "tukey")$se
amh_err_b <- mcse(mcmc.samp.b, method = "tukey")$se

df <- rbind(c(amh_mean_a, amh_mean_b), c(amh_err_a, amh_err_b))
xtable(df, digits = c(1,7,7))

#(c) Construct a variable-at-a-time Metropolis-Hastings (VMH) algorithm. You need to provide the 
# same level of detail here as you did for the previous algorithm.

proposal1=function(v1){
  return(rexp(1, rate = 1/v1))
}
proposal2=function(v2){
  return(rexp(1, rate = 1/v2))
}

target1=function(a,b){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(n*a*log(b) - b*sum.y - (1/6)*b^2)
    #return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y - (1/6)*a^2 - (1/6)*b^2)
    #return(log(prod(dgamma(y,a,b))*dnorm(b,0,3)))
  }
}


target2=function(a,b){
  if(alpha<0){
    return(-Inf)
  }else if(beta<0){
    return(-Inf)
  }
  else{
    return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - (1/6)*a^2)
    #return(n*a*log(b) - n*lgamma(a) + a*log.sum.y - b*sum.y - (1/6)*a^2- (1/6)*b^2)
    #return(log(prod(dgamma(y,a,b))*dnorm(a,0,3)))
  }
}


#Setting intial values and initializing variables
vmh.mcmc.samp.a=rep(NA,n.samp)
vmh.mcmc.samp.b=rep(NA,n.samp)
alpha=a.hat
beta=b.hat
vmh.mcmc.samp.a[1]=a.hat
vmh.mcmc.samp.b[1]=b.hat

#tuning parameters are just the means of the distribution

tic()
for(i in 1:n.samp){
  if(i %% 100 == 0) cat("Starting iteration", i, "\n")
  #VMH- beta
  beta.star=proposal1(beta)
  log.num=target1(alpha,beta.star)
  log.denom=target1(alpha,beta)
  r=exp(log.num-log.denom)
  if(runif(1)<r){
    beta=beta.star
  }
  
  #VMH- alpha
  alpha.star=proposal2(alpha)
  log.num=target2(alpha.star,beta)
  log.denom=target2(alpha,beta)
  r=exp(log.num-log.denom)
  if(runif(1)<r){
    alpha=alpha.star
  }
  
 vmh.mcmc.samp.a[i]=alpha
 vmh.mcmc.samp.b[i]=beta
}


vmh_time <- toc()
vmh_time <- vmh_time$toc-vmh_time$tic

## Diagnostics
# Acceptance rate
mean(!duplicated(vmh.mcmc.samp.a))
mean(!duplicated(vmh.mcmc.samp.b))

plot(vmh.mcmc.samp.a, type = "l", main ="Trace plots of alpha with burnin")
plot(vmh.mcmc.samp.b, type = "l", main ="Trace plots of beta with burnin")

burnin <- 2000
vmh.mcmc.samp.a <- vmh.mcmc.samp.a[-(1:burnin)]
vmh.mcmc.samp.b <- vmh.mcmc.samp.b[-(1:burnin)]

plot(vmh.mcmc.samp.a, type = "l", main = "Trace plots of alpha without burnin")
plot(vmh.mcmc.samp.b, type = "l", main = "Trace plots of beta without burnin")

acf(vmh.mcmc.samp.a, main = "ACF of alpha without burnin")
acf(vmh.mcmc.samp.b, main = "ACF of beta without burnin")

vmh_ess_a <- effectiveSize(vmh.mcmc.samp.a)
vmh_ess_b <- effectiveSize(vmh.mcmc.samp.b)

vmh_ess_a_t <- effectiveSize(vmh.mcmc.samp.a)/vmh_time
vmh_ess_b_t <- effectiveSize(vmh.mcmc.samp.b)/vmh_time

vmh_mean_a <- mean(vmh.mcmc.samp.a)
vmh_mean_b <- mean(vmh.mcmc.samp.b)
vmh_err_a <- mcse(vmh.mcmc.samp.a, method = "tukey")$se
vmh_err_b <- mcse(vmh.mcmc.samp.b, method = "tukey")$se

df <- rbind(c(vmh_mean_a, vmh_mean_b), c(vmh_err_a, vmh_err_b))
xtable(df, digits = c(1,7,7))

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
xtable(full_df, digits = c(7,7,7,7,7))

#(e) Compare the efficiency of the four algorithms {importance sampling, AMH, VMH, and the
# Laplace approximation from the previous homework} using the methodology discussed in class, 
# e.g. effective sample size, effective samples per second etc. The Monte Carlo algorithms are 
# easier to compare than the comparison between them and the Laplace approximation; you may 
# have to think carefully about this. Be sure to clearly explain why your approaches for comparing
# the algorithms are reasonable.
df <- cbind(c(imp_ess_a, imp_ess_b, imp_ess_a_t, imp_ess_b_t),
            c(amh_ess_a, amh_ess_b, amh_ess_a_t, amh_ess_b_t),
            c(vmh_ess_a, vmh_ess_b, vmh_ess_a_t, vmh_ess_b_t))
colnames(df) <- c("imp", "amh", "vmh")
xtable(df, digits=c(3,3,3,3))

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
proposal=function(x1){
  return(rexp(1, rate = 1/x1))
}

target=function(theta, z){
  if(theta<0){
    return(-Inf)
  }else{
    return( -100*log(theta) -(1/theta)*sum(bulb_dat) -(1/theta)*sum(z) -(theta/100))
  }
}

n.samp <- 10000
n <- length(bulb_dat)
mcmc.samp.aux=rep(NA,n.samp)
mcmc.samp.z <- matrix(ncol=25,nrow=n.samp)
mcmc.samp.z <- as.data.frame(mcmc.samp.z)
theta=100
mcmc.samp.aux[1]=theta
mcmc.samp.z[1,] <- r2exp(n=25, rate = theta, shift = 184)

tic()
for( i in 2:n.samp){
  if(i %% 1000 == 0) cat("Starting iteration", i, "\n")
  theta <- mcmc.samp.aux[i-1]
  x.star=proposal(theta)
  z <- mcmc.samp.z[i-1,]
  
  log.num=target(x.star,z)
  log.denom=target(theta,z)
  r=exp(log.num-log.denom)
  alpha <- min(1, r)
  
  if(runif(1)<alpha){
    theta=x.star
  }
  mcmc.samp.aux[i] <- theta
  
  z_val <- r2exp(n=25, rate = theta, shift = 184)
  mcmc.samp.z[i,] <- z_val
}
w_time <- toc()
w_time <- w_time$toc-w_time$tic

## Diagnostics
# Acceptance rate
mean(!duplicated(mcmc.samp.aux))

plot(mcmc.samp.aux, type = "l", main ="Trace plots of theta with burnin, with auxiliary")

burnin <- 500
mcmc.samp.aux <- mcmc.samp.aux[-(1:burnin)]

plot(mcmc.samp.aux, type = "l", main = "Trace plots of theta without burnin, with auxiliary")
acf(mcmc.samp.aux, main = "ACF of theta without burnin, with auxiliary")

ess_aux <- effectiveSize(mcmc.samp.aux)
es_t_aux <- effectiveSize(mcmc.samp.aux)/w_time

mean_aux <- mean(mcmc.samp.aux)
err_aux <- mcse(mcmc.samp.aux, method = "tukey")$se

df_aux <- as.data.frame(c(mean_aux, err_aux, ess_aux, es_t_aux))
xtable(df_aux, digits = c(7))

###############################################################################3
#b) Construct a different MCMC algorithm, this time without using auxiiliary variables/data 
# augmentation. Again, provide details.

# This time, I will approximate the posterior directly, without data augmentation.
 
proposal=function(x1){
  return(rexp(1, rate = 1/x1))
}

target=function(theta){
  if(theta<0){
    return(-Inf)
  }else{
    #return(-75*log(theta) -(1/theta)*sum(bulb_dat) - 75*log(1-exp(-(184/theta))) -(theta/100))
    return(-75*log(theta) -(1/theta)*sum(bulb_dat)-25*(184/theta)  - (theta/100))
  }
}


#Setting intial values and initializing variables
n <- length(bulb_dat)
mcmc.samp=rep(NA,n.samp)
theta=100
mcmc.samp[1]=theta

tic()
for(i in 2:n.samp){
  if(i %% 1000 == 0) cat("Starting iteration", i, "\n")
  theta <- mcmc.samp[i-1]
  x.star=proposal(theta)
  
  log.num=target(x.star)
  log.denom=target(theta)
  r=exp(log.num-log.denom)
  alpha <- min(1, r)
  
  if(runif(1)<alpha){
    theta=x.star
  }
  mcmc.samp[i]=theta
}
wo_time <- toc()
wo_time <- wo_time$toc-wo_time$tic

## Diagnostics
# Acceptance rate
mean(!duplicated(mcmc.samp))

plot(mcmc.samp, type = "l", main ="Trace plots of theta with burnin")

burnin <- 500
mcmc.samp <- mcmc.samp[-(1:burnin)]
mcmc.samp.2 <- mcmc.samp

plot(mcmc.samp, type = "l", main = "Trace plots of theta without burnin")
acf(mcmc.samp, main = "ACF of theta without burnin")

ess <- effectiveSize(mcmc.samp)
es_t <- effectiveSize(mcmc.samp)/wo_time

mean <- mean(mcmc.samp)
err <- mcse(mcmc.samp, method = "tukey")$se

df <- as.data.frame(c(mean, err, ess, es_t))
xtable(df, digits = c(7))

#(c) Overlay posterior density plot approximations for the two algorithms. Provide a table 
# that shows the posterior mean approximations for theta along with MCMC standard errors.
plot(density(mcmc_w_aux), "Density with auxiliarly variables")
lines(density(mcmc_wo_aux))

df <- cbind(c(mean_aux, err_aux, ess_aux, es_t_aux),c(mean, err, ess, es_t))

length(mcmc.samp.aux)
length(mcmc.samp)


#Sample data
dat <- data.frame(dens = c(mcmc.samp.aux, mcmc.samp),
                  lines = rep(c("with aux", "without aux"), each = length(mcmc.samp)))
#Plot.
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+labs(title ="Posterior Densities for two algorithms")


#(d) For the auxiliary variable method, plot the approximate posterior pdf for one of the 
# "missing" lightbulbs, then overlay the approximate posterior pdfs for a lightbulb made 
# by the company. You should notice that they are different. Report the posterior mean 
# estimates for each of them. 
plot(density(mcmc.samp.z[,1]), "missing lightbulb")

#just a lightbulb made by the company has a posterior (y|theta) of exponential, instead of shifted exponential
post_pdf_comp <- rexp(n=nrow(mcmc.samp.z),rate = 1/mean)

mean(mcmc.samp.z[,1])
mean(post_pdf_comp)

#Sample data
dat <- data.frame(dens = c(mcmc.samp.z[,1], post_pdf_comp),
                  lines = rep(c("missing lightbulb", "normal lightbulb"), each = length(post_pdf_comp)))
#Plot.
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+labs(title ="Missing vs Normal lightbulb")


#(e) Looking at just the ability to approximate the posterior per iteration of the algorithm, 
# which of the two MCMC algorithms is more efficient? Now accounting for computing costs, 
# which of the two MCMC algorithms is more efficient? Which algorithm would you recommend?
df <- cbind(c(ess_aux, es_t_aux),c(ess, es_t))
xtable(df, digits = c(7,7,7))

#(f) This course is focused on computing but it is worth noting some basics about inference. 
# Compare your results above with what would happen to inference if you ignored the missing 
# data by overlaying the density plots.

proposal=function(x1){
  return(rexp(1, rate = 1/x1))
}

target=function(theta){
  if(theta<0){
    return(-Inf)
  }else{
    return(-n*log(theta)-(1/theta)*(sum(bulb_dat)) +(theta/100))
    #return(dgamma(theta, shape = n+1, scale = (sum(bulb_dat)+ 1/100)^(-1), log = T))
  }
}


#Setting intial values and initializing variables
n <- length(bulb_dat)
mcmc.samp=rep(NA,n.samp)
theta=100
mcmc.samp[1]=theta

tic()
for(i in 2:n.samp){
  if(i %% 10000 == 0) cat("Starting iteration", i, "\n")
  theta <- mcmc.samp[i-1]
  x.star=proposal(theta)
  
  log.num=target(x.star)
  log.denom=target(theta)
  r=exp(log.num-log.denom)
  alpha <- min(1, r)
  
  if(runif(1)<alpha){
    theta=x.star
  }
  mcmc.samp[i]=theta
}
wo_time <- toc()
wo_time <- wo_time$toc-wo_time$tic

## Diagnostics
# Acceptance rate
mean(!duplicated(mcmc.samp))

plot(mcmc.samp, type = "l", main ="Trace plots of theta with burnin")

burnin <- 500
mcmc.samp <- mcmc.samp[-(1:burnin)]

plot(mcmc.samp, type = "l", main = "Trace plots of theta without burnin")
acf(mcmc.samp, main = "ACF of theta without burnin")

ess <- effectiveSize(mcmc.samp)
es_t <- effectiveSize(mcmc.samp)/wo_time

mean <- mean(mcmc.samp)
err <- mcse(mcmc.samp, method = "tukey")$se

df <- as.data.frame(c(mean, err, ess, es_t))
xtable(df, digits = c(7))


#Sample data
dat <- data.frame(dens = c(mcmc.samp.aux, mcmc.samp, mcmc.samp.2),
                  lines = rep(c("with aux", "ignoring aux", "without aux"), each = length(mcmc.samp)))
#Plot.
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)+labs(title ="Posterior Densities for three algorithms")
