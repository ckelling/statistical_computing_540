###
### Claire Kelling
### STAT 540
### Homework #1
###
### Created 2/21/2018 for the second assigment, due 3/8/2018
### The exercises mainly focus on ??
### 

library(MASS)
library(mvtnorm)
library(tictoc)
library(optimr)
library(ggplot2)
library(xtable)
library(tictoc)
library(RSpectra)

###
### Laplace Approximation
###
#clear workspace
rm(list=ls())

#load data, make sure to change to READ FROM WEBSITE *******
#http://personal.psu.edu/muh10/540/Rcode/hw1.dat
prob_5_dat <- read.table("C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/homework_1/data/prob_5_dat.txt", quote="\"", comment.char="")
y <- c((t(as.matrix(prob_5_dat))))
n <- length(y)

set.seed(123)

denom <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  denom <- log_post
}

num_b <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  num <- log_post + log(param["beta"])# + log(param["beta"])
}

num_a <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  num <- log_post + log(param["alpha"])# + log(param["beta"])
}

y <- rep(2,200)
#give a set of initial values
initial_values <- c(alpha = 3, beta = 0.5)
n_samp <- 10

tic()
exp_est_a <- rep(NA, n_samp)
exp_est_b <- rep(NA, n_samp)

for(i in 1:n_samp){
  #fnscale is -1 so that it maximizes the log posterior likelihood
  opt_fit_den <- optim(initial_values, denom, control = list(fnscale = -1), y = y, hessian = TRUE, method = c("L-BFGS-B"))
  opt_fit_num_a <- optim(initial_values, num_a, control = list(fnscale = -1), y = y, hessian = TRUE, method = c("L-BFGS-B"))
  opt_fit_num_b <- optim(initial_values, num_b, control = list(fnscale = -1), y = y, hessian = TRUE, method = c("L-BFGS-B"))
  
  
  full_den <- exp(opt_fit_den$value)*(det(solve(opt_fit_den$hessian))^(-1/2))
  full_num_a <- exp(opt_fit_num_a$value)*(det(solve(opt_fit_num_a$hessian))^(-1/2))
  full_num_b <- exp(opt_fit_num_b$value)*(det(solve(opt_fit_num_b$hessian))^(-1/2))
  
  exp_est_a[i] <- full_num_a/full_den
  exp_est_b[i] <- full_num_b/full_den
}
lap_time <- toc()
lap_time <- lap_time$toc-lap_time$tic
lap_mean_a <- mean(exp_est_a)
lap_err_a <- sd(exp_est_a)/sqrt(n_samp)
lap_mean_b <- mean(exp_est_b)
lap_err_b <- sd(exp_est_b)/sqrt(n_samp)


####
#### Problem 1
####

## (Return to the previous homework problem.) Let X1,..,Xn iid from Gamma(a,b) distribution. 
## Consider Bayesian inference for a,b with prior distributions N(0; 3) for both. 
## The data for this problem are available here: http://personal.psu.edu/muh10/540/Rcode/hw1.dat

#(a) Use importance sampling to approximate the expectations. Provide pseudocode, including relevant
# distributions. Is the vaability of your Monte Carlo approximation guaranteed to be finite? 
# Explain your answer. Provide Monte Criarlo standard errors for your estimates.



#(b) Construct an all-at-once Metropolis-Hastings (AMH) algorithm to approximate the posterior 
# expectation of a; b. Provide pseudocode for your algorithm. This should include any distributions 
# you had to derive. How did you determine the Monte Carlo approximations were good enough? 
# That is, discuss stopping criteria (how you determined chain length), the number of starting 
# values you tried, how you obtained initial values etc.


#(c) Construct a variable-at-a-time Metropolis-Hastings (VMH) algorithm. You need to provide the 
# same level of detail here as you did for the previous algorithm.


#(d) Provide a table with all approximations along with any error approximations, as well as the 
# computational time taken by the algorithms. (To make this easy to read, the rows of your table 
# should correspond to algorithms.)
lap_col <- c(lap_mean_a, lap_mean_b, lap_err_a, lap_err_b, lap_time)

full_df <- cbind(lap_col)
rownames(full_df) <- c("a_est", "b_est", "a_err", "b_err", "approx_time")


#(e) Compare the efficiency of the four algorithms {importance sampling, AMH, VMH, and the
# Laplace approximation from the previous homework} using the methodology discussed in class, 
# e.g. effective sample size, effective samples per second etc. The Monte Carlo algorithms are 
# easier to compare than the comparison between them and the Laplace approximation; you may 
# have to think carefully about this. Be sure to clearly explain why your approaches for comparing
# the algorithms are reasonable.


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
# (Assume that the remaining bulbs were still working at time tau.) Let the prior for tau be
# Gamma(1,100) (with parameterization so it has mean 100 and variance 100^2).

#load data
#http://personal.psu.edu/muh10/540/hwdir/bulbs.dat

#(a) Using auxiliary variables for the missing lightbulb data, construct an MCMC algorithm 
# to approximate the posterior distribution of tau. Provide the same level of detail about 
# your MCMC algorithm as you did in the previous problem.



#b) Construct a different MCMC algorithm, this time without using auxiiliary variables/data 
# augmentation. Again, provide details.


#(c) Overlay posterior density plot approximations for the two algorithms. Provide a table 
# that shows the posterior mean approximations for tau along with MCMC standard errors.
plot(density(df[,1]), "Density with auxiliarly variables")
lines(density(df[,2]), "Density without auxiliarly variables")

mean_aux <- mean(df[,1])
mean_wo_aux <- mean(df[,2])
se_aux <- sd(df[,1])/sqrt(n_samp)
se_wo_aux <- sd(df[,2])/sqrt(n_samp)

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
