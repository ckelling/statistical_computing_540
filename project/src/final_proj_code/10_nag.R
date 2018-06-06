###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/18/18 for final project code on comparison of Adam optimization (final algorithm)
### Here, we code SGD with momentum
### 
### Unlike for the algorithms, here we just create the storage structure for simplicity.
### 
library(data.table)
library(ggplot2)
library(tictoc)
library(xtable)


obj_fun = function(theta){
  x.b = as.vector(X%*%theta)
  return( -sum(y*x.b - log( 1 + exp(x.b))))  
  
}

Score.nag <- function(theta,batch.size, index){
  loop_sum = matrix(0,ncol=1,nrow=2)
  
  for(i in index){
    exp.i = exp(crossprod(X[i,],theta))
    const.i_num = as.numeric(exp.i)
    const.i_den = as.numeric((1+exp.i))
    
    loop_sum = loop_sum + X[i,] * (y[i] - const.i_num / const.i_den )
  }
  
  return(-loop_sum)
}




#
#Store multiple iterations of the algorithm:
#
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
    new_theta <- theta_prior - step*beta1*mean.curr
    Score.eval = Score.nag(new_theta,1, stoch_iter)
    
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
      Score.eval <- Score.nag(new_theta,1, stoch_iter)
      
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



##Load data and try algorithm
#load data
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")

maxit <- 500
eps2 <- 1e-6
beta1 <- 0.9
step <- 100
data <- prob1
init <- c(10,10)
data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]



niter <- 100
nag_store <- list()
nag_store$df <-NULL
nag_store$theta.final  <- NULL
nag_store$theta.1.hist  <- NULL
nag_store$theta.2.hist  <- NULL

#will need to delete this row at the end
nag_store$df <- rbind(nag_store$df, c(0,0,0,0,0,0))
i <- 0

while((nrow(nag_store$df) < niter)){
  i <- i + 1
  if(nrow(nag_store$df) %% 10 == 0) print(paste("we are currently at iteration", nrow(nag_store$df), " out of ", niter))
  
  #store time
  tic()
  nag_log_out <- nag_opt_st(data, step, eps2, beta1, init, maxit)
  nag_time <- toc()
  nag_time <- nag_time$toc-nag_time$tic
  
  #check for convergence
  th1 <- nag_log_out$theta.final[1]
  th2 <- nag_log_out$theta.final[2]
  b_iter <- nag_log_out$backtrack.iter
  
  #b_iter being greater than 20,000 means there is some problem with the starting/tuning 
  if((b_iter < 20000)){
    
    b_iter <- nag_log_out$backtrack.iter
    iter <- nag_log_out$iter
    conv_iter <- nag_log_out$conv_iter
    b_iter_c <- nag_log_out$b_iter_c[1]
    c_iter <- conv_iter[1]
    tot_iter <- sum(c_iter, b_iter_c)
    
    if(is.null(c_iter)){
      c_iter <- maxit
      b_iter_c <- b_iter
      print(paste("algo didn't converge at iteration ", i, " out of ", niter))
    }
    
    #store data on each iteration
    nag_store$df <- rbind(nag_store$df, c(nag_time, iter, b_iter, c_iter,b_iter_c,tot_iter))
    
    print(dim(nag_store$df))
    
    #numerical underflow problems, so will need to just include the converged number after a certain point
    filler_th1 <- rep(th1, 501-length(nag_log_out$theta.hist[,1]))
    filler_th2 <- rep(th2, 501-length(nag_log_out$theta.hist[,2]))
    
    #store data that will be dataframes
    nag_store$theta.final  <- rbind(nag_store$theta.final, nag_log_out$theta.final)
    nag_store$theta.1.hist  <- cbind(nag_store$theta.1.hist, c(nag_log_out$theta.hist[,1],filler_th1))
    nag_store$theta.2.hist  <- cbind(nag_store$theta.2.hist, c(nag_log_out$theta.hist[,2],filler_th2))
  }
}
nag_store$f_conv <- i-nrow(nag_store$df)
nag_store$perc_conv <- (i-nrow(nag_store$df))/i

nag_store$df <- nag_store$df[-1,]

colnames(nag_store$df) <- c("time", "iter", "b_iter", "c_iter","b_iter_c","tot_iter")

nag_theta_av <- NULL
nag_theta_av$th_1_av <- rowMeans(nag_store$theta.1.hist)
nag_theta_av$th_2_av <- rowMeans(nag_store$theta.2.hist)
nag_theta_av <- as.data.frame(nag_theta_av)
rownames(nag_theta_av) <- NULL
colnames(nag_theta_av) <- c("th1", "th2")

#plotting the averaged iteration
theta_av_df <- cbind(c(rep("nag", 2*nrow(nag_theta_av))),
                  c(rep("theta_1", nrow(nag_theta_av)),rep("theta_2", nrow(nag_theta_av))),
                  c(nag_theta_av$th1,nag_theta_av$th2 ))
theta_av_df <- as.data.frame(theta_av_df)
colnames(theta_av_df) <- c("algo","coeff", "est")
theta_av_df$ind <- c(1:nrow(nag_theta_av),1:nrow(nag_theta_av))
theta_av_df$est <- as.numeric(as.character(theta_av_df$est))

ggplot(data=theta_av_df, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line()+
  geom_point()+labs(title = "Convergence Comparion")+
  facet_wrap(~ coeff, ncol = 3)

nag_df <- theta_av_df 

#save(nag_store, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nag_sim.Rdata")
#save(nag_df, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nag_df.Rdata")


#load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nag_sim.Rdata")
