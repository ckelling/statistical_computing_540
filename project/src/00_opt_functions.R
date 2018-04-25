obj_fun = function(theta){
  x.b = as.vector(X%*%theta)
  return( -sum(y*x.b - log( 1 + exp(x.b))))  
  
}

Score.SGD <- function(theta,batch.size, index){
  batch.size <- 3000
  loop_sum = matrix(0,ncol=1,nrow=length(theta))
  
  for(i in index){
    #i <- 300
    exp.i = exp(crossprod(X[i,],theta))
    const.i_num = as.numeric(exp.i)
    const.i_den = as.numeric((1+exp.i))
    
    loop_sum = loop_sum + X[i,] * (y[i] - const.i_num / const.i_den )
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
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
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
    #print(t)
    
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
      Score.eval <- Score.SGD(theta_prior,1, stoch_iter)
      
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

nadam_st <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  
  #initialization of output format
  nadam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10; diff_stop <- T; step_stop <- T #arbitrarily high number and false
  
  #now, we will create the updates according to the algorithm
  #momentum scheduler
  i <- 1:(maxit+1)
  beta1.t.vec <- beta1*(1-0.5*0.96^(i/250))
  
  #while((diff_theta > eps2) & (t < maxit) & diff_stop & step_stop){
  while((diff_theta > eps2) & (t < maxit) & diff_stop & step_stop){
    #keep the prior theta to assess for convergence
    theta_prior <- theta_vec
    theta_vec <- as.vector(theta_vec)
    theta_prior <- as.vector(theta_prior)
    
    #update time step
    t <- t+1
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    # update theta
    #Get gradients w.r.t. stochastic objective at timestep t
    grad <- Score.eval
    g_hat <- (grad)/ (1-prod(beta1.t.vec[1:t]))
    
    #current estimate of theta
    theta <- theta_prior
    
    #Update biased first moment estimate
    mean <- beta1 * mean.curr + (1 - beta1) * grad
    #compute unbiased estimate, with moment scheduler
    mhat <- mean/(1-prod(beta1.t.vec[1:(t+1)]))
    
    #Update biased second raw moment estimate
    var <- beta2 * var.curr + (1 - beta2) * (grad^2)
    #compute unbiased estimate
    var_hat <- (var)/(1-beta2^t)
    
    #compute new update
    m_bar <- (1-beta1.t.vec[t])*g_hat + beta1.t.vec[t+1]*mhat
    
    #Compute bias-corrected first and second moment estimate
    #Update parameters
    theta_vec <- theta - step * m_bar / (sqrt(var_hat) + eps1)
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    #ll.new - ll.prev
    
    ## Backtracking (halve step size)- this happens once we get closer
    if(ll.new - ll.prev >= 0){
      backtrack = T
      s = step
      #bt_theta <- NULL
    }else{
      backtrack = F
    }
    while(backtrack){
      s = s / 2
      
      #find the stochastic row that we will use
      stoch_iter <- sample(1:nrow(data), 1)
      
      #need for calculating gradient for logistic regression
      Score.eval = Score.SGD(theta_prior,1, stoch_iter)
      
      # update theta
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval
      g_hat <- (grad)/ (1-prod(beta1.t.vec[1:t]))
      
      #Update biased first moment estimate
      mean <- beta1 * mean.curr + (1 - beta1) * grad
      #compute unbiased estimate, with moment scheduler
      mhat <- mean/(1-prod(beta1.t.vec[1:(t+1)]))
      
      #Update biased second raw moment estimate
      var <- beta2 * var.curr + (1 - beta2) * (grad^2)
      #compute unbiased estimate
      var_hat <- (var)/(1-beta2^t)
      
      #compute new update
      m_bar <- (1-beta1.t.vec[t])*g_hat + beta1.t.vec[t+1]*mhat
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta_vec <- theta - s * m_bar / (sqrt(var_hat) + eps1)
      
      ll.new = obj_fun(theta_vec)
      
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1
      #print(paste(backtrack.iter, "**********"))
      
      if(backtrack.iter >20000){
        #print(paste("did not converge"))
        break
      }
      
      #theta.hist <- rbind(theta.hist, c(theta_vec))
    }
    
    #theta.hist <- rbind(theta.hist, bt_theta)
    mean.curr <- mean
    var.curr <- var
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    #Create output vector
    nadam.hist   <- rbind(nadam.hist, c(t, lr_t, mean, var))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$nadam.hist   <- nadam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}


adam_st <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  
  
  #initialization of output format
  adam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
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
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    # update theta
    #Get gradients w.r.t. stochastic objective at timestep t
    grad <- Score.eval
    
    #current estimate of theta
    theta <- theta_prior
    
    #Update biased first moment estimate
    mean <- beta1 * mean.curr + (1 - beta1) * grad
    #Update biased second raw moment estimate
    var <- beta2 * var.curr + (1 - beta2) * (grad^2)
    
    # #compute lr_t
    # num <- 1 - (beta2^t)
    # denom <- 1 - (beta1^t)
    # lr_t <- step * sqrt(num)/denom
    # 
    # #Compute bias-corrected first and second moment estimate
    # #Update parameters
    # e_hat <- eps1*sqrt(1-beta2^t)
    # theta_vec <- theta - lr_t * mean / (sqrt(var)) #+ e_hat)
    
    #Now, instead of lr_t, we will compute bias_corrected versions of mean and variance, 
    #  and use these as the next estimates (this is not what is shown in the paper)
    m_hat <- mean/(1 - beta1^t)
    v_hat <- var/(1 - beta2^t)
    
    #Compute bias-corrected first and second moment estimate
    #Update parameters
    theta_vec <- theta_prior - step * m_hat / (sqrt(v_hat)+ eps1)
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    #ll.new - ll.prev
    
    ## Backtracking (halve step size)- this happens once we get closer
    if(ll.new - ll.prev >= 0){
      backtrack = T
      s = step
      #bt_theta <- NULL
    }else{
      backtrack = F
    }
    while(backtrack){
      s = s / 2
      
      #find the stochastic row that we will use
      stoch_iter <- sample(1:nrow(data), 1)
      
      #need for calculating gradient for logistic regression
      Score.eval = Score.SGD(theta_prior,1, stoch_iter)
      
      # update theta
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval
      
      #current estimate of theta
      theta <- theta_prior
      
      #Update biased first moment estimate
      mean <- beta1 * mean.curr + (1 - beta1) * grad
      #Update biased second raw moment estimate
      var <- beta2 * var.curr + (1 - beta2) * (grad^2)
      
      #Now, instead of lr_t, we will compute bias_corrected versions of mean and variance, 
      #  and use these as the next estimates (this is not what is shown in the paper)
      m_hat <- mean/(1 - beta1^t)
      v_hat <- var/(1 - beta2^t)
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta_vec <- theta - s * m_hat / (sqrt(v_hat) + eps1)
      
      ll.new = obj_fun(theta_vec)
      
      ll.new - ll.prev
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1
      #print(paste(backtrack.iter, "**********"))
      
      if(backtrack.iter >20000){
        #print(paste("did not converge"))
        break
      }
      
      #theta.hist <- rbind(theta.hist, c(theta_vec))
    }
    
    #theta.hist <- rbind(theta.hist, bt_theta)
    mean.curr <- mean
    var.curr <- var
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    
    #Create output vector
    adam.hist   <- rbind(adam.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$adam.hist   <- adam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}


nag_opt_st <- function(data, step, eps2, beta1, init, maxit){
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  theta_vec <- init
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  mean.curr <- 0
  
  #initialization of output format
  nag.hist <- data.frame(t)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
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
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    new_theta <- theta_prior - step*beta1*mean.curr
    Score.eval = Score.SGD(new_theta,1, stoch_iter)
    
    #momentum
    mean <- beta1*mean.curr + Score.eval
    
    # update theta
    theta_vec <- theta_prior - step * mean
    
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
      new_theta <- theta_prior - s*beta1*mean.curr
      Score.eval <- Score.SGD(new_theta,1, stoch_iter)
      
      #momentum
      mean <- beta1*mean.curr + Score.eval
      
      theta_vec <- theta_prior - s * mean  # update theta
      theta_vec <- as.vector(theta_vec)
      
      #check to see if continuing backtrack
      ll.new = obj_fun(theta_vec)
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1 
      
      #if(backtrack.iter %% 10000 == 0) print(backtrack.iter)
      #print(backtrack.iter)
      if(backtrack.iter > 20000) break
      
    }
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    
    #Create output vector
    nag.hist   <- rbind(nag.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    mean.curr <- mean
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$nag.hist   <- nag.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}


adam_st_nobc <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  
  
  #initialization of output format
  adam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
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
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    # update theta
    #Get gradients w.r.t. stochastic objective at timestep t
    grad <- Score.eval
    
    #current estimate of theta
    theta <- theta_prior
    
    #Update biased first moment estimate
    mean <- beta1 * mean.curr + (1 - beta1) * grad
    #Update biased second raw moment estimate
    var <- beta2 * var.curr + (1 - beta2) * (grad^2)
    
    # #compute lr_t
    # num <- 1 - (beta2^t)
    # denom <- 1 - (beta1^t)
    # lr_t <- step * sqrt(num)/denom
    # 
    # #Compute bias-corrected first and second moment estimate
    # #Update parameters
    # e_hat <- eps1*sqrt(1-beta2^t)
    # theta_vec <- theta - lr_t * mean / (sqrt(var)) #+ e_hat)
    
    #Now, instead of lr_t, we will compute bias_corrected versions of mean and variance, 
    #  and use these as the next estimates (this is not what is shown in the paper)
    #m_hat <- mean/(1 - beta1^t)
    #v_hat <- var/(1 - beta2^t)
    
    #Compute bias-corrected first and second moment estimate
    #Update parameters
    theta_vec <- theta_prior - step * mean / (sqrt(var)+ eps1)
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    #ll.new - ll.prev
    
    ## Backtracking (halve step size)- this happens once we get closer
    if(ll.new - ll.prev >= 0){
      backtrack = T
      s = step
      #bt_theta <- NULL
    }else{
      backtrack = F
    }
    while(backtrack){
      s = s / 2
      
      #find the stochastic row that we will use
      stoch_iter <- sample(1:nrow(data), 1)
      
      #need for calculating gradient for logistic regression
      Score.eval = Score.SGD(theta_prior,1, stoch_iter)
      
      # update theta
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval
      
      #current estimate of theta
      theta <- theta_prior
      
      #Update biased first moment estimate
      mean <- beta1 * mean.curr + (1 - beta1) * grad
      #Update biased second raw moment estimate
      var <- beta2 * var.curr + (1 - beta2) * (grad^2)
      
      #Now, instead of lr_t, we will compute bias_corrected versions of mean and variance, 
      #  and use these as the next estimates (this is not what is shown in the paper)
      #m_hat <- mean/(1 - beta1^t)
      #v_hat <- var/(1 - beta2^t)
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta_vec <- theta - s * mean / (sqrt(var) + eps1)
      
      ll.new = obj_fun(theta_vec)
      
      ll.new - ll.prev
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1
      #print(paste(backtrack.iter, "**********"))
      
      if(backtrack.iter >20000){
        #print(paste("did not converge"))
        break
      }
      
      #theta.hist <- rbind(theta.hist, c(theta_vec))
    }
    
    #theta.hist <- rbind(theta.hist, bt_theta)
    mean.curr <- mean
    var.curr <- var
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    
    #Create output vector
    adam.hist   <- rbind(adam.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$adam.hist   <- adam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}


nadam_st_nobc <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  
  #initialization of output format
  nadam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10; diff_stop <- T; step_stop <- T #arbitrarily high number and false
  
  #now, we will create the updates according to the algorithm
  #momentum scheduler
  i <- 1:(maxit+1)
  beta1.t.vec <- beta1*(1-0.5*0.96^(i/250))
  
  #while((diff_theta > eps2) & (t < maxit) & diff_stop & step_stop){
  while((diff_theta > eps2) & (t < maxit) & diff_stop & step_stop){
    #keep the prior theta to assess for convergence
    theta_prior <- theta_vec
    theta_vec <- as.vector(theta_vec)
    theta_prior <- as.vector(theta_prior)
    
    #update time step
    t <- t+1
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    # update theta
    #Get gradients w.r.t. stochastic objective at timestep t
    grad <- Score.eval
    g_hat <- (grad)/ (1-prod(beta1.t.vec[1:t]))
    
    #current estimate of theta
    theta <- theta_prior
    
    #Update biased first moment estimate
    mean <- beta1 * mean.curr + (1 - beta1) * grad
    #compute unbiased estimate, with moment scheduler
    #mhat <- mean/(1-prod(beta1.t.vec[1:(t+1)]))
    
    #Update biased second raw moment estimate
    var <- beta2 * var.curr + (1 - beta2) * (grad^2)
    #compute unbiased estimate
    #var_hat <- (var)/(1-beta2^t)
    
    #compute new update
    m_bar <- (1-beta1.t.vec[t])*g_hat + beta1.t.vec[t+1]*mean
    
    #Compute bias-corrected first and second moment estimate
    #Update parameters
    theta_vec <- theta - step * m_bar / (sqrt(var) + eps1)
    
    ll.prev = obj_fun(theta_prior)
    ll.new = obj_fun(theta_vec)
    
    #ll.new - ll.prev
    
    ## Backtracking (halve step size)- this happens once we get closer
    if(ll.new - ll.prev >= 0){
      backtrack = T
      s = step
      #bt_theta <- NULL
    }else{
      backtrack = F
    }
    while(backtrack){
      s = s / 2
      
      #find the stochastic row that we will use
      stoch_iter <- sample(1:nrow(data), 1)
      
      #need for calculating gradient for logistic regression
      Score.eval = Score.SGD(theta_prior,1, stoch_iter)
      
      # update theta
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- Score.eval
      g_hat <- (grad)/ (1-prod(beta1.t.vec[1:t]))
      
      #Update biased first moment estimate
      mean <- beta1 * mean.curr + (1 - beta1) * grad
      #compute unbiased estimate, with moment scheduler
      #mhat <- mean/(1-prod(beta1.t.vec[1:(t+1)]))
      
      #Update biased second raw moment estimate
      var <- beta2 * var.curr + (1 - beta2) * (grad^2)
      #compute unbiased estimate
      #var_hat <- (var)/(1-beta2^t)
      
      #compute new update
      m_bar <- (1-beta1.t.vec[t])*g_hat + beta1.t.vec[t+1]*mean
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta_vec <- theta - s * m_bar / (sqrt(var) + eps1)
      
      ll.new = obj_fun(theta_vec)
      
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1
      #print(paste(backtrack.iter, "**********"))
      
      if(backtrack.iter >20000){
        #print(paste("did not converge"))
        break
      }
      
      #theta.hist <- rbind(theta.hist, c(theta_vec))
    }
    
    #theta.hist <- rbind(theta.hist, bt_theta)
    mean.curr <- mean
    var.curr <- var
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    #Create output vector
    nadam.hist   <- rbind(nadam.hist, c(t, lr_t, mean, var))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$nadam.hist   <- nadam.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}


sgdm_opt_st <- function(data, step, eps2, beta1, init, maxit){
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  theta_vec <- init
  backtrack.iter <- 0
  conv_iter <- NULL
  b_iter_c <- NULL
  mean.curr <- 0
  
  #initialization of output format
  sgdm.hist <- data.frame(t)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data <- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- cbind(rep(1, nrow(data)), X) #also include an intercept
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
    if(t %% 1000 == 0) print(t)
    
    #find the stochastic row that we will use
    stoch_iter <- sample(1:nrow(data), 1)
    
    #need for calculating gradient for logistic regression
    Score.eval = Score.SGD(theta_prior,1, stoch_iter)
    
    #momentum
    mean <- beta1*mean.curr + Score.eval
    
    # update theta
    theta_vec <- theta_prior - step * mean
    
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
      Score.eval <- Score.SGD(theta_prior,1, stoch_iter)
      
      #momentum
      mean <- beta1*mean.curr + Score.eval
      
      theta_vec <- theta_prior - s * mean  # update theta
      theta_vec <- as.vector(theta_vec)
      
      #check to see if continuing backtrack
      ll.new = obj_fun(theta_vec)
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1 
      
      #if(backtrack.iter %% 10000 == 0) print(backtrack.iter)
      #print(backtrack.iter)
      if(backtrack.iter > 20000) break
      
    }
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    diff_theta <- max(abs(theta_prior-theta_vec)) > eps2
    
    if(diff_stop & ll_stop & diff_theta){
      c <- "do nothing"
    }else{
      conv_iter <- c(conv_iter,t)
      b_iter_c <- c(b_iter_c, backtrack.iter)
    }
    
    
    #Create output vector
    sgdm.hist   <- rbind(sgdm.hist, c(t))
    theta.hist <- rbind(theta.hist, c(theta_vec))
    
    mean.curr <- mean
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$conv_iter  <- conv_iter
  output$b_iter_c  <- b_iter_c
  output$iter        <- t - 1
  output$sgdm.hist   <- sgdm.hist
  output$theta.hist  <- theta.hist
  output$backtrack.iter <- backtrack.iter
  
  return(output)
}

