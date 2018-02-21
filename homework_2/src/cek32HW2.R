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
### Problem 1
###
# Let X1,...,Xn iid from Gamma(alpha; beta) distribution. Consider Bayesian
# inference for alpha; beta with prior distributions a~ N(0; 3), beta~ N(0; 3). Use a
# Laplace approximation (as discussed in class) to approximate the posterior
# expectation of alpha, beta.

#clear workspace
rm(list=ls())

#load data
prob_5_dat <- read.table("C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/homework_1/data/prob_5_dat.txt", quote="\"", comment.char="")
y <- c((t(as.matrix(prob_5_dat))))
n <- length(y)

set.seed(123)

# Find unnormalized log posterior of the model
#    given parameter vector pparam and vector of data points y



denom <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  denom <- log_post
}



num <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  num <- log_post + log(param["beta"])# + log(param["beta"])
}


#give a set of initial values
initial_values <- c(alpha = 3, beta = 0.5)

tic()
#fnscale is -1 so that it maximizes the log posterior likelihood
opt_fit_den <- optim(initial_values, denom, control = list(fnscale = -1), y = y, hessian = TRUE, method = c("L-BFGS-B"))
opt_fit_num <- optim(initial_values, num, control = list(fnscale = -1), y = y, hessian = TRUE, method = c("L-BFGS-B"))

full_den <- exp(opt_fit_den$value)*(det(solve(opt_fit_den$hessian))^(-1/2))
full_num <- exp(opt_fit_num$value)*(det(solve(opt_fit_num$hessian))^(-1/2))

exp_est <- full_num/full_den
toc()


##
## What if y has length 1,000? 2,000? Computation time?
##
prob5_t <- NULL
n <- c((1:10)*1000, (2:10)*10000, (2:10)*100000, (2:3)*1000000)
for (i in n){
  y <- runif(n = i, min = 1, max = 16)
  print(paste(i, "********************"))
  tic()
  #fnscale is -1 so that it maximizes the log posterior likelihood
  opt_fit_den <- optim(initial_values, denom, control = list(fnscale = -1), y = y, hessian = TRUE)
  print(paste("done with denominator"))
  opt_fit_num <- optim(initial_values, num, control = list(fnscale = -1), y = y, hessian = TRUE)
  
  full_den <- exp(opt_fit_den$value)*(det(solve(opt_fit_den$hessian))^(-1/2))
  full_num <- exp(opt_fit_num$value)*(det(solve(opt_fit_num$hessian))^(-1/2))
  
  exp_est <- full_num/full_den
  time <- toc()
  time <- time$toc-time$tic
  prob5_t <- c(prob5_t, time)
}

prob_5 <- data.frame(n, prob5_t)
colnames(prob_5) <- c("n", "system_time")


ggplot(prob_5) + geom_line(aes(x=n, y = system_time), size = 1.2)+
  labs(title = "Computational Cost scaling with n", y = "wall time", x = "length of y")
