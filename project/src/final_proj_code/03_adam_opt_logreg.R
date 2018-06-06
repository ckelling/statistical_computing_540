###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/20/18 for final project code on Adam optimization
### Here, we code the adam optimizer by hand, and assess for possible difficulties.
### This implementation is for the problem we did for Homework 3 (NR, GD, and SGD) for comparison
### In other words, this is trying to find the Gamma parameters
### 

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

adam_opt <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean <- 0
  var <- 0
  theta_vec <- init
  lr_t <- 0
  
  #initialization of output format
  adam.hist <- data.frame(t, lr_t, mean, var)
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
    for(k in 1:length(theta_vec)){
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval[k]
      
      #current estimate of theta
      theta <- theta_vec[k]
      
      #Update biased first moment estimate
      mean <- beta1 * mean + (1 - beta1) * grad
      #Update biased second raw moment estimate
      var <- beta2 * var + (1 - beta2) * (grad^2)
      
      #compute lr_t
      num <- 1 - beta1^t
      denom <- 1 - beta2^t
      lr_t <- step * sqrt(num)/denom
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta <- theta - lr_t * mean / (sqrt(var) + eps1)
      theta_vec[k] <- theta
    }
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    
    #Create output vector
    adam.hist   <- rbind(adam.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$iter        <- t - 1
  output$adam.hist   <- adam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}

#test the adam algorithm
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 500
lr <- 0.001
data <- prob1

adam_out <- adam_opt(data,lr, beta1, beta2, eps1, eps2, init, maxit)
theta_hist <- adam_out$theta.hist

##
## Now, I would like to apply Nesterov-accelerated adaptive moment estimation
##

nadam_opt <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean <- 0
  var <- 0
  theta_vec <- init
  
  #initialization of output format
  nadam.hist <- data.frame(t, mean, var)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data<- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10; diff_stop <- T; step_stop <- T #arbitrarily high number and false
  
  #momentum scheduler
  i <- 1:maxit
  beta1.t.vec <- beta1*(1-0.5*0.96^(i/250))
  
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
    for(k in 1:length(theta_vec)){
      #current estimate of theta
      theta <- theta_vec[k]
      
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval[k]
      g_hat <- (grad)/ (1-prod(beta1.t.vec[1:t]))
      
      #Update biased first moment estimate
      mean <- beta1 * mean + (1 - beta1) * grad
      #compute unbiased estimate, with moment scheduler
      mhat <- mean/(1-prod(beta1.t.vec[1:(t+1)]))
      
      #Update biased second raw moment estimate
      var <- beta2 * var + (1 - beta2) * (grad^2)
      #compute unbiased estimate
      var_hat <- (var)/(1-beta2^t)
      
      #compute new update
      m_bar <- (1-beta1.t.vec[t])*g_hat + beta1.t.vec[t+1]*mhat
      
      #Update parameters
      theta <- theta - lr * m_bar / (sqrt(var_hat) + eps1)
      theta_vec[k] <- theta
    }
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    
    #Create output vector
    nadam.hist   <- rbind(nadam.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$iter        <- t - 1
  output$nadam.hist   <- nadam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}

#Running Nadam algorithm
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 700
lr <- 0.001

nadam_output <- nadam_opt(data,lr, beta1, beta2, eps1, eps2, init, maxit)
theta_hist <- nadam_output$theta.hist


###
### Measurements for Adam
###

tic()
adam_out <- adam_opt(data,0.1, beta1, beta2, eps1, eps2, init, maxit)
adam_time <- toc()
adam_time <- adam_time$toc-adam_time$tic

adam_out$iter

theta_df <- cbind(c(rep("theta_1", nrow(adam_out$theta.hist)),rep("theta_2", nrow(adam_out$theta.hist))),
                  c(adam_out$theta.hist[,1],adam_out$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(adam_out$theta.hist),1:nrow(adam_out$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  #geom_line()+
  geom_point()+labs(title = "Adam Algorithm")

adam_out$theta.final

###
### Measurements for Nadam
###

tic()
nadam_out <- nadam_opt(data,0.1, beta1, beta2, eps1, eps2, init, maxit)
nadam_time <- toc()
nadam_time <- nadam_time$toc-nadam_time$tic

nadam_out$iter

theta_df <- cbind(c(rep("theta_1", nrow(nadam_out$theta.hist)),rep("theta_2", nrow(nadam_out$theta.hist))),
                  c(nadam_out$theta.hist[,1],nadam_out$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(nadam_out$theta.hist),1:nrow(nadam_out$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  #geom_line()+
  geom_point()+labs(title = "Nadam Algorithm")


