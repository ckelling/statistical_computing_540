###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/20/18 for final project code on nadam optimization
### Here, we code the nadam optimizer by hand, and assess for possible difficulties.
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

nadam_bc <- function(data, step, beta1, beta2, eps1, eps2, init, maxit){  
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean.curr <- c(0,0)
  var.curr <- c(0,0)
  theta_vec <- init
  lr_t <- 0
  backtrack.iter <- 0
  
  #initialization of output format
  nadam.hist <- data.frame(t, lr_t, mean.curr, var.curr)
  theta.hist <- theta_vec
  
  #split data into x and y variables
  data<- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
  X <- as.matrix(X)
  y <- data[,ncol(data)]
  diff_theta <- 10; diff_stop <- T; step_stop <- T #arbitrarily high number and false
  
  #now, we will create the updates according to the algorithm
  #momentum scheduler
  i <- 1:maxit
  beta1.t.vec <- beta1*(1-0.5*0.96^(i/250))
  
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
    nadam.hist   <- rbind(nadam.hist, c(t, lr_t, mean, var))
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

#test the nadam algorithm
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 10000
step <- 1#step size of 1 helps speed
init <- c(10,10)
data <- prob1

data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]

###
### Measurements for nadam
###

#tic()
nadam_bc_out <- nadam_bc(data, step, beta1, beta2, eps1, eps2, init, maxit)
theta_hist <- nadam_bc_out$theta.hist
#nadam_time <- toc()
#nadam_time <- nadam_time$toc-nadam_time$tic

nadam_bc_out$iter

theta_df <- cbind(c(rep("theta_1", nrow(nadam_bc_out$theta.hist)),rep("theta_2", nrow(nadam_bc_out$theta.hist))),
                  c(nadam_bc_out$theta.hist[,1],nadam_bc_out$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(nadam_bc_out$theta.hist),1:nrow(nadam_bc_out$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  #geom_line()+
  geom_point()+labs(title = "Nadam Algorithm")

nadam_bc_out$backtrack.iter
nadam_bc_out$theta.final


#save(nadam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_output.Rdata")
#save(nadam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_output2.Rdata")


###
### Now I will create a version of the algorithm that will make it easier to store information
###
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
  data<- as.data.frame(data)
  X <- data[,1:(ncol(data)-1)]
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

##Load data and try algorithm for storing
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 500
step <- 1#step size of 1 helps speed
init <- c(10,10)
data <- prob1

data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]


niter <- 100
nadam_store <- list()
nadam_store$df <-NULL
nadam_store$theta.final  <- NULL
nadam_store$theta.1.hist  <- NULL
nadam_store$theta.2.hist  <- NULL

#will need to delete this row at the end
nadam_store$df <- rbind(nadam_store$df, c(0,0,0,0,0,0))
j <- 0

while((nrow(nadam_store$df) < niter)){
  j <- j + 1
  if(nrow(nadam_store$df) %% 10 == 0) print(paste("we are currently at iteration", nrow(nadam_store$df), " out of ", niter))
  
  #store time
  tic()
  nadam_log_out <- nadam_st(data, step, beta1, beta2, eps1, eps2, init, maxit)
  #nadam_log_out$theta.final
  #b_iter <- nadam_log_out$backtrack.iter
  #b_iter
  adam_time <- toc()
  adam_time <- adam_time$toc-adam_time$tic
  
  #check for convergence
  th1 <- nadam_log_out$theta.final[1]
  th2 <- nadam_log_out$theta.final[2]
  
  if((th1 < -0.6 & th1 >-1) & (th2 > 1.7 & th2 <2)){
    
    b_iter <- nadam_log_out$backtrack.iter
    iter <- nadam_log_out$iter
    conv_iter <- nadam_log_out$conv_iter
    b_iter_c <- nadam_log_out$b_iter_c[1]
    c_iter <- conv_iter[1]
    tot_iter <- sum(c_iter, b_iter_c)
    
    if(is.null(c_iter)){
      c_iter <- maxit
      b_iter_c <- b_iter
      print(paste("algo didn't converge at iteration ", i, " out of ", niter))
    }
    
    #store data on each iteration
    nadam_store$df <- rbind(nadam_store$df, c(adam_time, iter, b_iter, c_iter,b_iter_c,tot_iter))
    
    print(dim(nadam_store$df))
    
    #numerical underflow problems, so will need to just include the converged number after a certain point
    filler_th1 <- rep(th1, 501-length(nadam_log_out$theta.hist[,1]))
    filler_th2 <- rep(th2, 501-length(nadam_log_out$theta.hist[,2]))
    #store data that will be dataframes
    nadam_store$theta.final  <- rbind(nadam_store$theta.final, nadam_log_out$theta.final)
    nadam_store$theta.1.hist  <- cbind(nadam_store$theta.1.hist, c(nadam_log_out$theta.hist[,1],filler_th1))
    nadam_store$theta.2.hist  <- cbind(nadam_store$theta.2.hist, c(nadam_log_out$theta.hist[,2],filler_th2))
  }
}
nadam_store$f_conv <- j-nrow(nadam_store$df)
nadam_store$perc_conv <- (j-nrow(nadam_store$df))/j

#save(nadam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_output.Rdata")
#save(nadam_bc_out, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_output2.Rdata")

nadam_store$df <- nadam_store$df[-1,]

colnames(nadam_store$df) <- c("adam_time", "iter", "b_iter", "c_iter","b_iter_c","tot_iter")

nadam_theta_av <- NULL
nadam_theta_av$th_1_av <- rowMeans(nadam_store$theta.1.hist)
nadam_theta_av$th_2_av <- rowMeans(nadam_store$theta.2.hist)
nadam_theta_av <- as.data.frame(nadam_theta_av)
rownames(nadam_theta_av) <- NULL
colnames(nadam_theta_av) <- c("th1", "th2")

#plotting the averaged iteration
theta_av_df <- cbind(c(rep("nadam", 2*nrow(nadam_theta_av))),
                     c(rep("theta_1", nrow(nadam_theta_av)),rep("theta_2", nrow(nadam_theta_av))),
                     c(nadam_theta_av$th1,nadam_theta_av$th2 ))
theta_av_df <- as.data.frame(theta_av_df)
colnames(theta_av_df) <- c("algo","coeff", "est")
theta_av_df$ind <- c(1:nrow(nadam_theta_av),1:nrow(nadam_theta_av))
theta_av_df$est <- as.numeric(as.character(theta_av_df$est))

ggplot(data=theta_av_df, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line()+
  geom_point()+labs(title = "Convergence Comparion")+
  facet_wrap(~ coeff, ncol = 3)