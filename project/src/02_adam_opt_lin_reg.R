###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/20/18 for final project code on Adam optimization
### Here, we code the adam optimizer by hand, and assess for possible difficulties.
### This implementation is for linear regression.
### 


### This package allows for parallelization.
# cran <- getOption("repos")
# cran["dmlc"] <- "https://apache-mxnet.s3-accelerate.dualstack.amazonaws.com/R/CRAN/"
# options(repos = cran)
# install.packages("mxnet")
# install.packages("https://github.com/jeremiedb/mxnet_winbin/raw/master/mxnet.zip", repos = NULL)

library(mxnet)
data(gradDescentRData)

# This is a dataset that collected experimentaly by Kennedy in 1954 to obtain the density value of
# CO2.  Parameter used in the experiment are temperature and pressure, which it can be a parameter
# to obtain compressibility factor value of CO2.

#first, try fitting linear model and these can be initial values for theta
## get z-factor data
dataSet <- gradDescentRData$CompressilbilityFactor

## split dataset into training and testing data
split_dat <- splitData(dataSet)
tr_dat <- split_dat$dataTrain
data <- tr_dat
colnames(data) <- c("temp", "pressure", "comp_fac")
ln_mod <- lm(comp_fac ~ temp + pressure, data = data)
init <- ln_mod$coefficients


#function for calculating the L2 norm
norm_vec <- function(x) sqrt(sum(x^2))


adam_opt <- function(data, lr, beta1, beta2, eps1, eps2, init, maxit){
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean <- 0
  var <- 0
  theta_vec <- init
  lr_t <- 0
  
  #initialization of output format
  adam.hist <- data.frame(t, lr_t, mean, var)
  theta.hist <- theta_vec
  
  #give data an intercept
  data <- cbind(1, data)
  
  #split data into x and y variables
  x_var <- data[,1:ncol(data)-1]
  x_var <- as.matrix(x_var)
  y_var <- data[,ncol(data)]
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
    error <- (t(x_var[stoch_iter,]) %*% (theta_vec)) - y_var[stoch_iter]
    
    for(k in 1:length(theta_vec)){
      #Get gradients w.r.t. stochastic objective at timestep t
      #k <- 2
      grad <- error * x_var[stoch_iter, k]
      
      #current estimate of theta
      theta <- theta_vec[k]
      
      #Update biased first moment estimate
      mean <- beta1 * mean + (1 - beta1) * grad
      #Update biased second raw moment estimate
      var <- beta2 * var + (1 - beta2) * (grad^2)
      
      #compute lr_t
      num <- 1 - beta1^t
      denom <- 1 - beta2^t
      lr_t <- lr * sqrt(num)/denom
      
      #Compute bias-corrected first and second moment estimate
      #Update parameters
      theta <- theta - lr_t * mean / (sqrt(var) + eps1)
      theta_vec[k] <- theta
    }
    
    #Create output vector
    adam.hist   <- rbind(adam.hist, c(t, lr_t, mean, var))
    theta.hist <- rbind(theta.hist, theta_vec)
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$iter        <- t - 1
  output$adam.hist   <- adam.hist
  output$theta.hist  <- theta.hist
  
  return(output)
  
}

#try it out!
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 0.0001 #convergence criteria
maxit <- 500
lr <- 0.001

adam_out <- adam_opt(tr_dat,lr, beta1, beta2, eps1, eps2, init, maxit)
theta_hist <- test$theta.hist

##
## Now, I would like to apply Nesterov-accelerated adaptive moment estimation
##

nadam_opt <- function(data, lr, beta1, beta2, eps1, eps2, init, maxit){
  
  #initializations of iteration counts and other variables, as in paper
  t <- 0
  mean <- 0
  var <- 0
  theta_vec <- init
  
  #initialization of output format
  nadam.hist <- data.frame(t, mean, var)
  theta.hist <- theta_vec
  
  #give data an intercept
  #data <- cbind(1, tr_dat) #need to change back to data
  data <- cbind(1, data)
  
  #split data into x and y variables
  x_var <- data[,1:ncol(data)-1]
  x_var <- as.matrix(x_var)
  y_var <- data[,ncol(data)]
  diff_theta <- 10
  
  #momentum scheduler
  i <- 1:maxit
  beta1.t.vec <- beta1*(1-0.5*0.96^(i/250))
  
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
    error <- (t(x_var[stoch_iter,]) %*% (theta_vec)) - y_var[stoch_iter]
    
    for(k in 1:length(theta_vec)){
      #for testing purposes
      #k<-2
      
      #current estimate of theta
      theta <- theta_vec[k]
      
      #Get gradients w.r.t. stochastic objective at timestep t
      grad <- error * x_var[stoch_iter, k]
      
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
    
    #Create output vector
    nadam.hist   <- rbind(adam.hist, c(t, mhat, var_hat))
    theta.hist <- rbind(theta.hist, theta_vec)
    
    diff_theta <- max(abs(theta_prior-theta_vec))
  }
  
  output <- list()
  output$theta.final  <- theta_vec
  output$iter        <- t - 1
  output$adam.hist   <- nadam.hist
  output$theta.hist  <- theta.hist
  
  return(output)
  
}

#testing it out!
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 0.0001 #convergence criteria
maxit <- 500
lr <- 0.001

test <- nadam_opt(tr_dat,lr, beta1, beta2, eps1, eps2, init, maxit)
theta_hist <- test$theta.hist


###
### Measurements for Adam
###

tic()
adam_out <- adam_opt(tr_dat,lr, beta1, beta2, eps1, eps2, init, maxit)
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

###
### Measurements for Nadam
###

tic()
nadam_out <- nadam_opt(tr_dat,lr, beta1, beta2, eps1, eps2, init, maxit)
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