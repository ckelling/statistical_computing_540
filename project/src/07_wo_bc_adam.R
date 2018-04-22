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
    
    loop_sum = loop_sum + X[i,] * (y[i] - const.i_num / const.i_den )
  }
  
  return(-loop_sum)
}

adam_bc <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  
  #initialization of output format
  adam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
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
    
    #Now, instead of lr_t, we will compute bias_corrected versions of mean and variance, 
    #  and use these as the next estimates (this is not what is shown in the paper)
    #mean <- mean/(1 - beta1^t)
    #var <- var/(1 - beta2^t)
    
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
      #mean <- mean/(1 - beta1^t)
      #var <- var/(1 - beta2^t)

      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta_vec <- theta - s * mean / (sqrt(var) + eps1)

      ll.new = obj_fun(theta_vec)

      ll.new - ll.prev
      backtrack = ll.new - ll.prev >= 0
      backtrack.iter = backtrack.iter + 1
      #print(paste(backtrack.iter, "**********"))

      if(backtrack.iter >40000){
        print(paste("did not converge"))
        break
      }

      #theta.hist <- rbind(theta.hist, c(theta_vec))
    }

    #theta.hist <- rbind(theta.hist, bt_theta)
    mean.curr <- mean
    var.curr <- var
    
    diff_stop = max( abs(theta_prior-theta_vec)/ abs(theta_prior)) > 1e-10  
    ll_stop   = abs( ll.new-ll.prev ) / abs(ll.prev)  > 1e-10 
    
    #Create output vector
    adam.hist   <- rbind(adam.hist, c(t, lr_t, mean, var))
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
maxit <- 10000
step <- 1#step size of 1 helps speed (can make step size smaller for more consistent convergence)
init <- c(10,10)
data <- prob1

data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]


###
### Measurements for Adam
###

#tic()
adam_bc_out <- adam_bc(data, step, beta1, beta2, eps1, eps2, init, maxit)
#theta_hist <- adam_bc_out$theta.hist
#adam_time <- toc()
#adam_time <- adam_time$toc-adam_time$tic

adam_bc_out$iter

theta_df <- cbind(c(rep("theta_1", nrow(adam_bc_out$theta.hist)),rep("theta_2", nrow(adam_bc_out$theta.hist))),
                  c(adam_bc_out$theta.hist[,1],adam_bc_out$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(adam_bc_out$theta.hist),1:nrow(adam_bc_out$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  #geom_line()+
  geom_point()+labs(title = "Adam Algorithm")

adam_bc_out$backtrack.iter
adam_bc_out$theta.final

#save(adam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/adam_output.Rdata")
#save(adam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/adam_output3.Rdata")


###
### Now I will create a version of the algorithm that will make it easier to store information
###
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

##Load data and try algorithm for storing
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 500
step <- 1#step size of 1 helps speed (can make step size smaller for more consistent convergence)
init <- c(10,10)
data <- prob1

data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]


niter <- 100
adam_store <- list()
adam_store$df <-NULL
adam_store$theta.final  <- NULL
adam_store$theta.1.hist  <- NULL
adam_store$theta.2.hist  <- NULL

#will need to delete this row at the end
adam_store$df <- rbind(adam_store$df, c(0,0,0,0,0,0))
i <- 0

while((nrow(adam_store$df) < niter)){
  i <- i + 1
  if(nrow(adam_store$df) %% 10 == 0) print(paste("we are currently at iteration", nrow(adam_store$df), " out of ", niter))
  
  #store time
  tic()
  adam_log_out <- adam_st(data, step, beta1, beta2, eps1, eps2, init, maxit)
  adam_log_out$theta.final
  b_iter <- adam_log_out$backtrack.iter
  b_iter
  adam_time <- toc()
  adam_time <- adam_time$toc-adam_time$tic
  
  #check for convergence
  th1 <- adam_log_out$theta.final[1]
  th2 <- adam_log_out$theta.final[2]
  
  if((th1 < -0.6 & th1 >-1) & (th2 > 1.7 & th2 <2) & (b_iter< 20000)){
    
    b_iter <- adam_log_out$backtrack.iter
    iter <- adam_log_out$iter
    conv_iter <- adam_log_out$conv_iter
    b_iter_c <- adam_log_out$b_iter_c[1]
    c_iter <- conv_iter[1]
    tot_iter <- sum(c_iter, b_iter_c)
    
    if(is.null(c_iter)){
      c_iter <- maxit
      b_iter_c <- b_iter
      print(paste("algo didn't converge at iteration ", i, " out of ", niter))
    }
    
    #store data on each iteration
    adam_store$df <- rbind(adam_store$df, c(adam_time, iter, b_iter, c_iter,b_iter_c,tot_iter))
    
    print(dim(adam_store$df))
    
    #numerical underflow problems, so will need to just include the converged number after a certain point
    filler_th1 <- rep(th1, 501-length(adam_log_out$theta.hist[,1]))
    filler_th2 <- rep(th2, 501-length(adam_log_out$theta.hist[,2]))
    
    #store data that will be dataframes
    adam_store$theta.final  <- rbind(adam_store$theta.final, adam_log_out$theta.final)
    adam_store$theta.1.hist  <- cbind(adam_store$theta.1.hist, adam_log_out$theta.hist[,1])
    adam_store$theta.2.hist  <- cbind(adam_store$theta.2.hist, adam_log_out$theta.hist[,2])
  }
}
adam_store$f_conv <- i-nrow(adam_store$df)
adam_store$perc_conv <- (i-nrow(adam_store$df))/i

###
### Save this data and plot it for just the algorithm
###

adam_store$df <- adam_store$df[-1,]

colnames(adam_store$df) <- c("adam_time", "iter", "b_iter", "c_iter","b_iter_c","tot_iter")

adam_theta_av <- NULL
adam_theta_av$th_1_av <- rowMeans(adam_store$theta.1.hist)
adam_theta_av$th_2_av <- rowMeans(adam_store$theta.2.hist)
adam_theta_av <- as.data.frame(adam_theta_av)
rownames(adam_theta_av) <- NULL
colnames(adam_theta_av) <- c("th1", "th2")

#plotting the averaged iteration
theta_av_df <- cbind(c(rep("adam", 2*nrow(adam_theta_av))),
                     c(rep("theta_1", nrow(adam_theta_av)),rep("theta_2", nrow(adam_theta_av))),
                     c(adam_theta_av$th1,adam_theta_av$th2 ))
theta_av_df <- as.data.frame(theta_av_df)
colnames(theta_av_df) <- c("algo","coeff", "est")
theta_av_df$ind <- c(1:nrow(adam_theta_av),1:nrow(adam_theta_av))
theta_av_df$est <- as.numeric(as.character(theta_av_df$est))

ggplot(data=theta_av_df, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line()+
  geom_point()+labs(title = "Convergence Comparison")+
  facet_wrap(~ coeff, ncol = 3)

#save(adam_store_nobc, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_adam_sim_nobc.Rdata")
