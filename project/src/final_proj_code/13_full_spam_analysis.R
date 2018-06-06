###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/15/18 for final project use of spam dataset
### 

library(data.table)
library(ggplot2)
library(tictoc)
library(xtable)

#load the data

#https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.DOCUMENTATION
#https://archive.ics.uci.edu/ml/datasets/spambase
#https://www.kaggle.com/c/cs-189-logistic-regression2/data
spam_dat <- fread("https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data")

#load all of the optimization functions
source(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/src/00_opt_functions.R")

#Below, is our final dataset for our report, with 6 covariates (including an intercept)
data <- spam_dat[,c(1,2,3,5,6,58)]

#test fitting the glm to see if there is a fit
mod <- glm(V58 ~.,family=binomial(link='logit'),data=data)
summary(mod)

#values for the algorithm
beta1 <- 0.9
beta1 <- 0.99
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-9 #convergence criteria
maxit <- 3000
step <- 1#step size of 1 helps speed

#structuring the data, and adding an intercept
data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- cbind(rep(1, nrow(data)), X) #also include an intercept
X <- as.matrix(X)
y <- data[,ncol(data)]

#intial values for the algorithm
init <- rep(1, ncol(X))

### We substitute these values for the algorithm in the following while loop
# sgd_spam <- sgd_st(data, 100, eps2, init, maxit)
# sgd_spam$conv_iter
# sgdm_spam <- sgdm_opt_st(data, step, eps2, beta1, init, maxit)
# nag_spam <- nag_opt_st(data, step, eps2, beta1, init, maxit)
# adam_spam <- adam_st(data, 0.1, beta1, beta2, eps1, eps2, init, maxit)
# nadam_spam <- nadam_st(data, step, beta1, beta2, eps1, eps2, init, maxit)
# adam_nobc_spam <- adam_st_nobc(data, step, beta1, beta2, eps1, eps2, init, maxit)
# nadam_nobc_spam <- nadam_st_nobc(data, step, beta1, beta2, eps1, eps2, init, maxit)


###
### Structure for saving information
### 
niter <- 10
store <- list()
store$df <-NULL
store$theta.final  <- NULL
store$theta.1.hist  <- NULL
store$theta.2.hist  <- NULL
store$obj_st <- NULL
step <- 0.1 #0.1
#step <- 10
#will need to delete this row at the end
store$df <- rbind(store$df, c(0,0,0,0,0,0))
i <- 0

#for sgd, we use step size of 10
#for sgdm, we use step size of 0.002
#for adam, use 0.1
while((nrow(store$df) < niter)){
  i <- i + 1
  if(nrow(store$df) %% 10 == 0) print(paste("we are currently at iteration", nrow(store$df), " out of ", niter))
  
  #store time
  tic()
  spam_out <- adam_st(data, 0.1, beta1, beta2, eps1, eps2, init, maxit)
  sgdm_time <- toc()
  sgdm_time <- sgdm_time$toc-sgdm_time$tic
  
  #check for convergence
  th1 <- spam_out$theta.final[1]
  th2 <- spam_out$theta.final[2]
  b_iter <- spam_out$backtrack.iter
  conv_iter <- spam_out$conv_iter
  c_iter <- conv_iter[1]
  iter <- spam_out$iter
  
  print(spam_out$theta.final[1:2])
  print(b_iter)
  print(spam_out$iter)
  
  #once again, check for convergence (sometimes the algorithm gets stuck in backtracking)
  ##b_iter being greater than 20,000 means there is some problem with the starting/tuning
  if((b_iter < 20000)){
    
    b_iter <- spam_out$backtrack.iter
    iter <- spam_out$iter
    conv_iter <- spam_out$conv_iter
    b_iter_c <- spam_out$b_iter_c[1]
    c_iter <- conv_iter[1]
    tot_iter <- sum(c_iter, b_iter_c)
    
    if(is.null(c_iter)){
      c_iter <- maxit
      b_iter_c <- b_iter
      print(paste("algo didn't converge at iteration ", i, " out of ", niter))
    }
    
    #store data on each iteration
    store$df <- rbind(store$df, c(sgdm_time, iter, b_iter, c_iter,b_iter_c,tot_iter))
    
    print(dim(store$df))
    
    obj_final <- obj_fun(spam_out$theta.final)
    
    #numerical underflow problems, so will need to just include the converged number after a certain point
    #after convergence criteria has been satisfied, we keep this value for the rest (so the time is accurate)
    filler_th1 <- rep(th1, maxit+1-length(spam_out$theta.hist[,1]))
    filler_th2 <- rep(th2, maxit+1-length(spam_out$theta.hist[,2]))
    filler_obj <- rep(obj_final, maxit+1-length(spam_out$theta.hist[,2]))
    
    #store data that will be dataframes
    store$theta.final  <- cbind(store$theta.final, spam_out$theta.final)
    store$theta.1.hist  <- cbind(store$theta.1.hist, c(spam_out$theta.hist[,1],filler_th1))
    store$theta.2.hist  <- cbind(store$theta.2.hist, c(spam_out$theta.hist[,2],filler_th2))
    
    obj_fun_st <- rep(NA, nrow(spam_out$theta.hist))
    #information on convergence of objective function
    for(i in 1:nrow(spam_out$theta.hist)){
      obj_fun_st[i] <- obj_fun(spam_out$theta.hist[i,])
    }
    store$obj_st <- cbind(store$obj_st, c(obj_fun_st,filler_obj))
  }
}

#keep track of how many iterations fail
store$f_conv <- i-nrow(store$df)
store$perc_conv <- (i-nrow(store$df))/i

store$df <- store$df[-1,]

colnames(store$df) <- c("time", "iter", "b_iter", "c_iter","b_iter_c","tot_iter")

theta_av <- NULL
theta_av$th_1_av <- rowMeans(store$theta.1.hist)
theta_av$th_2_av <- rowMeans(store$theta.2.hist)
theta_av <- as.data.frame(theta_av)
rownames(theta_av) <- NULL
colnames(theta_av) <- c("th1", "th2")

#plotting the averaged iteration
theta_av_df <- cbind(c(rep("adam", 2*nrow(theta_av))),
                     c(rep("theta_1", nrow(theta_av)),rep("theta_2", nrow(theta_av))),
                     c(theta_av$th1,theta_av$th2 ))
theta_av_df <- as.data.frame(theta_av_df)
colnames(theta_av_df) <- c("algo","coeff", "est")
theta_av_df$ind <- c(1:nrow(theta_av),1:nrow(theta_av))
theta_av_df$est <- as.numeric(as.character(theta_av_df$est))

ggplot(data=theta_av_df, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line()+
  geom_point()+labs(title = "Convergence Comparion")+
  facet_wrap(~ coeff, ncol = 3)

rowMeans(store$theta.final)[1:2]
# -1.1413010  0.7394682 for sgd, could be used as a benchmark
# -1.1423518  0.7160199 second time around

#rename these as the file that we will save below
a_store2 <- store
a_storage_df2 <- theta_av_df

#sgd
#save(s_store2,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/try2/spam_full_sgd_sim.Rdata")
#save(s_storage_df2, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/try2/spam_sgd_df.Rdata")
#sgd with momentum
#save(sm_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_sgdm_sim.Rdata")
#save(sm_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_sgdm_df.Rdata")
#nag
#save(nag_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nag_sim.Rdata")
#save(nag_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nag_df.Rdata")
#adam
#save(a_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_sim2.Rdata")
#save(a_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_adam_df.Rdata")
#nadam
#save(n_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_sim.Rdata")
#save(n_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nadam_df.Rdata")
#adam, no bias correction
#save(abc_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_adam_nobc_sim.Rdata")
#save(abc_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_adam_df_nobc.Rdata")
#nadam, no bias correction
#save(nbc_store,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_nobc_sim.Rdata")
#save(nbc_storage_df,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nadam_df_nobc.Rdata")

