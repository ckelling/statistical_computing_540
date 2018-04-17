###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/28/18 for final project code on comparison of Adam optimization
### Here, we code SGD for comparison between Adam, Nadam, and SGD
### This implementation is for the problem we did for Homework 3 (NR, GD, and SGD) for comparison
### 
library(data.table)
library(ggplot2)
library(tictoc)
library(xtable)


obj_fun = function(theta){
  x.b = as.vector(X%*%theta)
  return( -sum(y*x.b - log( 1 + exp(x.b))))  
  
}

Score.SGD <- function(theta,batch.size, index){
  loop_sum = matrix(0,ncol=1,nrow=2)
  
  for(i in index){
    exp.i = exp(crossprod(X[i,],theta))
    const.i_num = as.numeric(exp.i)
    const.i_den = as.numeric((1+exp.i))
    
    loop_sum = loop_sum + X[i,] * (Y[i] - const.i_num / const.i_den )
  }
  
  return(-loop_sum)
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
  data<- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10; diff_stop <- T; step_stop <- T #arbitrarily high number and false
  
  #now, we will create the updates according to the algorithm
  while((diff_theta > eps2) & (t < maxit) & diff_stop & step_stop){
    
    #keep the prior theta to assess for convergence
    theta_prior <- theta_vec
    theta_vec <- as.vector(theta_vec)
    theta_prior <- as.vector(theta_prior)
    
    #update time step
    t <- t+1
    print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    # update theta
    theta_vec <- theta_prior - step * Score.eval
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    ### Backtracking (halve step size)
    if(ll.new - ll.prev >= 0){
      backtrack = T 
      s = step
    }else{
      backtrack = F
    }
    while(backtrack){
      s = s / 2
      
      #find the stochastic row that we will use
      stoch_iter <- sample(1:nrow(data), 1)
      Score.eval = Score.SGD(theta_prior,1, stoch_iter)
      
      theta_vec <- theta_prior - s * Score.eval  # update theta
      theta_vec <- as.vector(theta_vec)
      
      #print(paste("problem here   2"))
      ll.new = obj_fun(theta_vec)
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1 
      
      #print(backtrack.iter)

    }
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    
    #Create output vector
    sgd.hist   <- rbind(sgd.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
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

maxit <- 500
eps2 <- 1e-6
step <- 100
data <- prob1
init <- c(10,10)


tic()
sgd_log_out <- sgd_opt(data, step, eps2, init, maxit)
SGD_time <- toc()
SGD_time <- SGD_time$toc-SGD_time$tic

sgd_log_out$iter

theta_hist <- sgd_log_out$theta.hist

theta_df <- cbind(c(rep("theta_1", nrow(sgd_log_out$theta.hist)),rep("theta_2", nrow(sgd_log_out$theta.hist))),
                  c(sgd_log_out$theta.hist[,1],sgd_log_out$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(sgd_log_out$theta.hist),1:nrow(sgd_log_out$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  geom_line()+
  geom_point()+labs(title = "Stochastic Gradient Descent")

sgd_log_out$theta.final
