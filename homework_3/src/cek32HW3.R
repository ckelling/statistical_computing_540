###
### Claire Kelling
### STAT 540
### Homework #3
###
### Created 3/27/2018 for the third assigment, due 4/5/2018
### The exercises mainly focus on optimization algorithms
### 

library(data.table)
library(ggplot2)
library(tictoc)
library(xtable)

###
### Problem 1
###

# Find the MLE for a logistic regression model.
# Write a Newton-Raphon algorithm to find the MLE (\hat{\beta_1}; \hat{\beta_2}). 
# The data for this problem are available here: http://personal.psu.edu/muh10/540/data/logReg.dat

#clear workspace
rm(list=ls())

#load data
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")

##########################################
## Newton-Raphson
##########################################
#Find the MLE for a logistic regression model.} For $i = 1,...,n$
#  $$ Y_i \sim  Ber(p_i(\theta)); where  \theta= (\beta_1, \beta_2)$$
#    $$p_i(\theta) = exp(\beta_1X_{1i} + \beta_2X_{2i})/(1 + exp(\beta_1X_{1i} + \beta_2X_{2i}))$$ 
# where each $Y_i$ is a binary response corresponding to predictors $X_{1i},X_{2i}$. 
# Write a Newton-Raphson algorithm to find the MLE $( \hat{\beta}_1, \hat{\beta}_2)$. 

set.seed(1)
#plot(prob1$x1, prob1$y, main="x1 vs y")
#plot(prob1$x2, prob1$y, main="x2 vs y")

# creating function for probability (first derivative)
p_funct <- function(X, beta){
  # compute vector p of probabilities for logistic regression with logit link
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

# creating function for log-likelihood
log_lik <- function(y, p){
  # bernoulli log likelihood function
  # input:   vectors: y = 0/1 outcome; p = probabilities
  # output: log-likelihood l, a scalar
  one_vec <- rep(1,nrow(prob1))
  l <- t(y) %*% log(p) + t(one_vec - y) %*% log(1 - p)
  return(l)
}


#Now, I will create a function for Newton-Raphson, that includes some of the stopping criteria
f_NR <- function(X, y, beta.1, eps1, eps2, eps3, maxit, step){
  #input:
  #    (X,y): data, covariates and response
  #    beta.init: starting values for regression
  #    eps1: convergence criterion for beta
  #    eps2: convergence criterion for log-likelihood
  #    maxit: number of maximum iterations
  
  #initial beta.2, largely irrelevant
  beta.2 <- rep(-Inf, length(beta.1))
  
  #record the initial difference in beta and difference in log-likelihood (irrelevant)
  #   beta
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  
  #   log-likelihood
  llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood, beta.1
  llike.2 <- log_lik(y, p_funct(X, beta.2)) # update loglikelihood, beta.2
  diff.like <- abs(llike.1 - llike.2) # diff
  if (is.nan(diff.like)){
    diff.like <- 1e9
  }
  
  if (is.nan(llike.1)){
    llike.1 <- 1e9
  }
  score.2 <- t(X) %*% (y - p_funct(X, beta.1)) 
  #making structure to keep the data to analyze when it is done
  i <- 1
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1,score.2[1],score.2[2], step.size =1)
  beta.hist <- matrix(beta.1, nrow = 1)
  
  #actual Newton-Raphson iterations:
  #   while we have not had convergence and are not at our maximum number of iterations....
  while((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2) &
        (abs(score.2[1]) > eps3) & (abs(score.2[2]) > eps3)){
    # update index
    i <- i + 1                       
    
    #define step
    #step <- c(1,1)
    
    # update beta
    beta.2 <- beta.1
    p.2 <-  p_funct(X, beta.2)
    
    # score function
    score.2 <- t(X) %*% (y - p.2)   
    
    # variance matrix
    v.2  <- diag(as.vector(p_funct(X, beta.2) * (1 - p_funct(X, beta.2))))
    
    # this increment version inverts the information matrix
    # Iinv.2 <- solve(t(X) %*% v.2 %*% X)  # Inverse information matrix
    # increm <- Iinv.2 %*% score.2     # increment, solve() is inverse
    # this increment version solves for (beta.2-beta.1) without inverting Information
    increm <- solve(t(X) %*% v.2 %*% X, score.2) # solve for increment

    #update beta
    beta.1 <- beta.2 + step * increm   # update beta
    
    #now, calculate new distances (to assess for convergence)
    diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
    llike.2 <- llike.1               # age likelihood value
    llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood
    #if (is.nan(llike.1)){
    #  llike.1 <- 1e10
    #}
    
    
    diff.like <- abs(llike.1 - llike.2) # diff
    
    # iteration history
    NR.hist   <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1, score.2[1], score.2[2], step))
    beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
  }
  
  # prepare output
  #       beta.MLE: after convergence, what is beta_MLE
  #       iter: how long it took to converge
  #       NR.hist: history matrix as above, all information on convergence differences
  #       beta.hist: history of betas through iterations
  output <- list()
  output$beta.MLE  <- beta.1
  output$iter      <- i - 1
  output$NR.hist   <- NR.hist
  output$beta.hist <- beta.hist
  
  
  v.1  <- diag(as.vector(p_funct(X, beta.1) * (1 - p_funct(X, beta.1))))
  Iinv.1 <- solve(t(X) %*% v.1 %*% X)  # Inverse information matrix
  output$beta.cov  <- Iinv.1
  
  return(output)
}

#Now, try running the algorithm
n <- nrow(prob1)
y <- prob1$y
x1 <- prob1$x1
x2 <- prob1$x2
# X matrix
X <- matrix(c(x1, x2), nrow = n)
colnames(X) <- c("x1", "x2")
r <- ncol(X)

# initial beta vector
beta.zero <- rep(0, r)

# Change starting values to GLM fit
glm_fit <- glm(y~.-1,data.frame(y,x1,x2),family="binomial") 
beta.1.glm <- c(glm_fit$coefficients[1], glm_fit$coefficients[2])

#off coefficients
beta.off <- c(2,2)

#vary the step size
step_half <- c(0.5, 0.5)
step_1  <- c(1,1)
step_vary <- c(1.5,0.5)

# fit betas using our Newton Raphson function
tic()
output_0 <- f_NR(X, y, beta.zero, 1e-6, 1e-6, 1e-6, 50, step_vary)
NR_time <- toc()
NR_time <- NR_time$toc-NR_time$tic

output_0

beta_df <- cbind(c(rep("beta_1", nrow(output_0$beta.hist)),rep("beta_2", nrow(output_0$beta.hist))),
                 c(output_0$beta.hist[,1],output_0$beta.hist[,2] ))
colnames(beta_df) <- c("coeff", "est")
beta_df <- as.data.frame(beta_df)
beta_df$ind <- c(1:nrow(output_0$beta.hist),1:nrow(output_0$beta.hist))
beta_df$est <- as.numeric(as.character(beta_df$est))

ggplot(data=beta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  geom_line()+
  geom_point()+labs(title = "Step Size (1.5,0.5)")



# (a) Provide details about how you set up the algorithm, including any
#  analytical quantities you derived. You should also include the stopping criteria, 
#  how you obtained initial values, how you tuned sn (step size/learning rate) and other 
#  relevant information, e.g. did you run the algorithm using multiple initial values?

# (b) Provide your pseudocode for the main algorithm. You can use notation defined in class 
#  (and specified above) in order to make the pseudocode succinct and easy to read. 

# (c) Provide the MLE as well as estimate of the standard error, along
#  with asymptotic 95% confidence intervals.
beta.est  <- output_0$beta.MLE
beta.se   <- sqrt(diag(output_0$beta.cov))
CI_l <- beta.est - 1.96*beta.se
CI_u <- beta.est + 1.96*beta.se

results <- data.frame(beta.est, beta.se, CI_l, CI_u)
results

xtable(results, digits = c(5,5,5,5,5))

# (d) Report the total computing time taken by your algorithm (e.g. you
#  could use system.time). If you have 2-3 different versions of the
#  algorithm, report timings for each version.
NR_time
NR_time1

# (e) Provide the computational cost of the algorithm (flops) per iteration
#  as a function of n. Also provide the number of iterations it took before
#  the stopping criteria was met. Of course, the number of iterations will
#  vary by initial conditions so you should provide at least 2-3 different
#  counts depending on where you started.

# (f) Summarize in 1-2 sentences the sensitivity of the algorithm to the
#  initial value, that is, does the algorithm work well regardless of where
#  you start or does it need to be in a neighborhood of a certain value?
#    (This is not required, but because this is a simple 2-D optimization
#     problem, you could also provide the 2-D log-likelihood surface to
#     obtain insights.)



###
### Problem 2
###

# Implement a gradient descent algorithm (work with the negative log like-lihood). 
# Answer (a)-(f) for this algorithm.

#Now, I will create a function for gradient descent, that includes some of the stopping criteria
f_gd <- function(X, y, beta.1, eps1, eps2, eps3, maxit, step){
  #input:
  #    (X,y): data, covariates and response
  #    beta.init: starting values for regression
  #    eps1: convergence criterion for beta
  #    eps2: convergence criterion for log-likelihood
  #    maxit: number of maximum iterations
  
  #initial beta.2, largely irrelevant
  beta.2 <- rep(-Inf, length(beta.1))
  
  #record the initial difference in beta and difference in log-likelihood (irrelevant)
  #   beta
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  
  #   log-likelihood
  llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood, beta.1
  llike.2 <- log_lik(y, p_funct(X, beta.2)) # update loglikelihood, beta.2
  diff.like <- abs(llike.1 - llike.2) # diff
  if (is.nan(diff.like)){
    diff.like <- 1e9
  }
  
  #making structure to keep the data to analyze when it is done
  score.2 <- t(X) %*% (y - p_funct(X, beta.1)) 
  i <- 1
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1,score.2[1],score.2[2], step)
  beta.hist <- matrix(beta.1, nrow = 1)
  
  #actual gradient descent iterations:
  #   while we have not had convergence and are not at our maximum number of iterations....
  while((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2) &
        (abs(score.2[1]) > eps3) & (abs(score.2[2]) > eps3)){
    # update index
    i <- i + 1                       
    
    # update beta
    beta.2 <- beta.1
    p.2 <-  p_funct(X, beta.2)
    
    # score function
    score.2 <- t(X) %*% (y - p.2)   

    # this increment version uses the score function (or the gradient of f)
    increm <-  -score.2/length(y) # solve for increment
    
    #update beta
    beta.1 <- beta.2 - step * increm   # update beta
    
    #now, calculate new distances (to assess for convergence)
    diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
    llike.2 <- llike.1               # likelihood value
    llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood
    diff.like <- abs(llike.1 - llike.2) # diff
    
    # iteration history
    NR.hist   <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1, score.2[1], score.2[2], step))
    beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
  }
  
  # prepare output
  #       beta.MLE: after convergence, what is beta_MLE
  #       iter: how long it took to converge
  #       NR.hist: history matrix as above, all information on convergence differences
  #       beta.hist: history of betas through iterations
  output <- list()
  output$beta.MLE  <- beta.1
  output$iter      <- i - 1
  output$NR.hist   <- NR.hist
  output$beta.hist <- beta.hist
  
  v.1  <- diag(as.vector(p_funct(X, beta.1) * (1 - p_funct(X, beta.1))))
  Iinv.1 <- solve(t(X) %*% v.1 %*% X)  # Inverse information matrix
  output$beta.cov  <- Iinv.1
  
  return(output)
}

#Now, try running the algorithm

# initial beta vector choices
beta.zero <- rep(0, r)
beta.1.glm <- c(glm_fit$coefficients[1], glm_fit$coefficients[2])
beta.off <- c(2,2)
beta.close <- c(-1,2)

#vary the step size
step_half <- c(0.5, 0.5)
step_1  <- c(1,1)
step_vary <- c(1.5,0.5)
step_ten <- c(10,10)
step_15 <- c(15,15)


# fit betas using our Gradient Descent function
tic()
output_0 <- f_gd(X, y, beta.zero, 1e-6, 1e-6, 1e-6, 700, step_ten)
gd_time <- toc()
gd_time <- gd_time$toc-gd_time$tic

#output_0
output_0$iter

beta_df <- cbind(c(rep("beta_1", nrow(output_0$beta.hist)),rep("beta_2", nrow(output_0$beta.hist))),
                 c(output_0$beta.hist[,1],output_0$beta.hist[,2] ))
colnames(beta_df) <- c("coeff", "est")
beta_df <- as.data.frame(beta_df)
beta_df$ind <- c(1:nrow(output_0$beta.hist),1:nrow(output_0$beta.hist))
beta_df$est <- as.numeric(as.character(beta_df$est))

ggplot(data=beta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  geom_line()+
  geom_point()+labs(title = "step (10,10), initial values (-1,2)")


# Estimates: 
beta.est  <- output_0$beta.MLE
beta.se   <- sqrt(diag(output_0$beta.cov))
CI_l <- beta.est - 1.96*beta.se
CI_u <- beta.est + 1.96*beta.se

results <- data.frame(beta.est, beta.se, CI_l, CI_u)
results

xtable(results, digits = c(5,5,5,5,5))




###
### Problem 3
###

# Now implement a stochastic gradient descent algorithm and repeat (a)-(f).

#Now, I will create a function for gradient descent, that includes some of the stopping criteria
f_sgd <- function(X, y, beta.1, eps1, eps2,eps3, maxit, step){
  #input:
  #    (X,y): data, covariates and response
  #    beta.init: starting values for regression
  #    eps1: convergence criterion for beta
  #    eps2: convergence criterion for log-likelihood
  #    maxit: number of maximum iterations
  
  #initial beta.2, largely irrelevant
  beta.2 <- rep(-Inf, length(beta.1))
  
  #record the initial difference in beta and difference in log-likelihood (irrelevant)
  #   beta
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  
  #   log-likelihood
  llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood, beta.1
  llike.2 <- log_lik(y, p_funct(X, beta.2)) # update loglikelihood, beta.2
  diff.like <- abs(llike.1 - llike.2) # diff
  if (is.nan(diff.like)){
    diff.like <- 1e9
  }
  
  #making structure to keep the data to analyze when it is done
  score.2 <- t(X) %*% (y - p_funct(X, beta.1)) 
  i <- 1
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1,score.2[1],score.2[2], step)
  beta.hist <- matrix(beta.1, nrow = 1)
  #j <- sample(seq(1,n), 1)
  
  #actual gradient descent iterations:
  #   while we have not had convergence and are not at our maximum number of iterations....
  while((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2) &
        (abs(score.2[1]) > eps3) & (abs(score.2[2]) > eps3)){
    # update index
    i <- i + 1                       
    
    # update beta
    beta.2 <- beta.1
    
    p.2 <- p_funct(X,beta.2)
    #p.2 <- exp(beta.2 %*% t(X[j,]))/(1+exp(beta.2 %*% t(X[j,])))
    
    score.1 <- score.2
    
    # score function
    j <- sample(seq(1,n), 1) 
    #j <- i %% nrow(X) + 2
    score.2 <- t(X[j,]) * (y[j] - p.2[j])
    
    # this increment version uses the score function (or the gradient of f)
    increm <-  -score.2 #/n
    
    # update beta
    beta.1 <- beta.2 - step * increm   # update beta
    
    #now, calculate new distances (to assess for convergence)
    diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
    llike.2 <- llike.1               # likelihood value from last iteration
    llike.1 <- log_lik(y, p_funct(X, beta.1)) # update loglikelihood
    diff.like <- abs(llike.1 - llike.2) # diff
    
    #step size schedule
    if(llike.1 > llike.2){
      step <- 0.5*step
    }
    
    
    # iteration history
    NR.hist   <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1,score.2[1],score.2[2], step))
    beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
  }
  
  # prepare output
  #       beta.MLE: after convergence, what is beta_MLE
  #       iter: how long it took to converge
  #       NR.hist: history matrix as above, all information on convergence differences
  #       beta.hist: history of betas through iterations
  output <- list()
  output$beta.MLE  <- beta.1
  output$iter      <- i - 1
  output$NR.hist   <- NR.hist
  output$beta.hist <- beta.hist
  
  
  v.1  <- diag(as.vector(p_funct(X, beta.1) * (1 - p_funct(X, beta.1))))
  Iinv.1 <- solve(t(X) %*% v.1 %*% X)  # Inverse information matrix
  output$beta.cov  <- Iinv.1
  
  return(output)
}

# initial beta vector choices
beta.zero <- rep(0, r)
beta.1.glm <- c(glm_fit$coefficients[1], glm_fit$coefficients[2])
beta.off <- c(2,2)
beta.cl <- c(-2,2)
beta.close <- c(-1,2)
beta.close.2 <- c(-0.75,1.75)

#vary the step size
step_half <- c(0.5, 0.5)
step_1  <- c(1,1)
step_vary <- c(1.5,0.5)
step_ten <- c(10,10)
step_15 <- c(15,15)
step_05 <- c(0.05, 0.05)
step_100 <- c(100,100)


# fit betas using our SGD function
tic()
output_0 <- f_sgd(X, y, c(0,0), 1e-6, 1e-6, 1e-6, 10000, c(10,10))
SGD_time <- toc()
SGD_time <- SGD_time$toc-SGD_time$tic

output_0$iter

beta_df <- cbind(c(rep("beta_1", nrow(output_0$beta.hist)),rep("beta_2", nrow(output_0$beta.hist))),
                 c(output_0$beta.hist[,1],output_0$beta.hist[,2] ))
colnames(beta_df) <- c("coeff", "est")
beta_df <- as.data.frame(beta_df)
beta_df$ind <- c(1:nrow(output_0$beta.hist),1:nrow(output_0$beta.hist))
beta_df$est <- as.numeric(as.character(beta_df$est))

ggplot(data=beta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  geom_line()+
  geom_point()+labs(title = "step (0.5,0.5), initial values (2,2)")


# Estimates: 
beta.est  <- output_0$beta.MLE
beta.se   <- sqrt(diag(output_0$beta.cov))
CI_l <- beta.est - 1.96*beta.se
CI_u <- beta.est + 1.96*beta.se

results <- data.frame(beta.est, beta.se, CI_l, CI_u)
results

xtable(results, digits = c(5,5,5,5,5))


###
### Problem 4
###

#  Compare the three algorithms above: Newton-Raphson, gradient descent,
# and stochastic gradient descent. Provide a 1-2 sentence summary of which
# algorithm you would recommend and why. Then provide a more detailed
# comparison, for example how stable the algorithms are with respect to
# initial values, how sensitive they are to choice of s_n, a comparison of the
# total computational cost, and the computational cost per iteration.


###
### Problem 5
###

# E-M algorithm: Return to a lightbulb lifetimes problem similar to the
# problem in a previous homework. Assume lightbulbs made by a par-
#   ticular company are independent and gamma distributed with param-
#   eters alpha; beta (parameterization is such that expected value is alpha*beta). Sup-
#   pose in an experiment m bulbs are switched on at the same time, but
# are only completely observed up to time tau . Let the lifetimes of these
# bulbs be A1;A2; : : : ;Am. However, since the bulbs are only observed
# till time tau, not all these lifetimes will be observed. Now suppose that
# at time tau, the experimenter observes the number of lightbulbs, W, still
# working at time tau , and the lifetimes of all lightbulbs that stopped work-
#   ing by tau . For convenience, denote these bulbs that stopped working by
# time tau as A1;A2; : : : ;Am????W. Hence, the missing information consists of
# the lifetimes of the lightbulbs still working at time tau , Am????W+1; : : : ;Am.
# For a particular experiment, let tau be 200 days and m = 300. The
# data on the lightbulb lifetimes for the bulbs that stopped working by
# tau are here: http://personal.psu.edu/muh10/540/data/bulbsHW3.dat
# Assume that the remaining bulbs were still working at time tau . Find the
# MLE for (alpha; beta) using the E-M algorithm. Repeat parts (a)-(f) from above
# for this algorithm as well.

times <- fread("http://personal.psu.edu/muh10/540/data/bulbsHW3.dat")
