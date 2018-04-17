###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/28/18 for final project code on comparison of Adam optimization
### Here, we code SGD for comparison between Adam, Nadam, and SGD
### This implementation is for the problem we did for Homework 3 (NR, GD, and SGD) for comparison
### 

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


sgd_opt <- function(data, step, eps2, init, maxit){
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  theta_vec <- init
  backtrack.iter <- 0
  
  #initialization of output format
  sgd.hist <- data.frame(t)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  X <- data[,1:ncol(data)-1]
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10
  
  #now, we will create the updates according to the algorithm
  while((diff_theta > eps2) & (t < maxit)){
    
    #keep the prior theta to assess for convergence
    theta_prior <- theta_vec
    
    #update time step
    t <- t+1
    print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for linear regression
    #theta_vec <- as.numeric(theta_vec$theta_vec)
    p.2 <- p_funct(X,beta.2)
    
    score.1 <- score.2
    
    # score function
    # now, try to have a bigger batch size
    score.2 <- t(X[stoch_iter,]) %*% (y[stoch_iter] - p.2[stoch_iter])
    
    # this increment version uses the score function (or the gradient of f)
    increm <-  -score.2 /length(j)
    
    # update beta
    theta_vec <- theta_prior - step * increm   # update theta
    
    ll.prev = log_lik(y, p_funct(X, theta_prior))
    ll.new = log_lik(y, p_funct(X, theta_vec))
    
    ### Backtracking (halve step size)
    if(ll.new - ll.prev >= 0){
      backtrack = T 
      s = step.size
    }else{
      backtrack = F
    }
    while(backtrack){
      step = step / 2
      
      p.2 <- p_funct(X,theta_prior)
      score.2 <- t(X[stoch_iter,]) %*% (y[stoch_iter] - p.2[stoch_iter])
      increm <-  -score.2 
      
      beta.new = beta.prev - step*increm
      ll.new = Objective_fun(beta.new)
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1 
    }
    
    #Create output vector
    sgd.hist   <- rbind(sgd.hist, c(t))
    theta.hist <- rbind(theta.hist, theta_vec)
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$iter        <- t - 1
  output$sgd.hist   <- sgd.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}



##Load data and try algorithm
#load data
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
