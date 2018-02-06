###
### Claire Kelling
### STAT 540
### Homework #1
###
### Created 1/29/2018 
### 

library(MASS)
library(mvtnorm)
library(tictoc)
library(optimr)
library(ggplot2)
library(xtable)
library(tictoc)
library(RSpectra)

###
### Problem 1
###

# Write computer code to simulate M = 100 draws from a multivariate normal
#     distribution with mean 0 (vector of 0s) and covariance matrix
#     as above for n = 1000.
#     Sigma_{ij} = exp(-|i-j|/phi),

#Creating 3 different algorithms

sim_alg1 <- function(n, mu, Sigma){ #cholesky factorization
  X <- matrix(rnorm(n*M), nrow = M, ncol=n)
  L <- chol(Sigma)
  mv_n <- X%*%L #+ matrix(mu, nrow = M, ncol =n),  mu is 0 in our case so omitted
  return(list(mv_n,L))
}

sim_alg2 <- function(n, mu, Sigma){ #spectral decomposition
  dim <- length(mu)
  ev <- eigen(Sigma)
  val <- ev$values
  vec <- ev$vectors
  
  R <- vec%*%diag(sqrt(val))%*%t(vec)
  X <- matrix(rnorm(n*M), nrow = M, ncol=n)
  
  mv_n2 <- X%*%R + matrix(mu, nrow = M, ncol =n)
  return(mv_n2)
}

sim_alg3 <- function(n, mu, Sigma){ #SVD
  dim <- length(mu)
  S_vd <- svd(Sigma)
  u <- S_vd$u
  d <- S_vd$d
  v <- S_vd$v
  
  R <- u%*%diag(sqrt(d))%*%t(v)
  X <- matrix(rnorm(n*M), nrow = M, ncol=n)
  
  mv_n3 <- X%*%R + matrix(mu, nrow = M, ncol =n)
  return(mv_n3)
}
#https://books.google.com/books?id=cR1jDAAAQBAJ&pg=PA71&dq=generating+multivariate+normal+in+r&hl=en&sa=X&ved=0ahUKEwiCsaOW24PZAhVIj1kKHdQxDAsQ6AEIJzAA#v=onepage&q&f=false


#Write your own code to also evaluate the multivariate normal pdf at each simulated value.
#   https://www.cs.cmu.edu/~epxing/Class/10701-08s/recitation/gaussian.pdf
#   mu is 0 in this case so I will take it out

#Calculate pdf for each X

mvn_pdf <- function(X,C){
  #C <- chol(Sigma)
  
  #as in class, r'r = (x-mu)Sigma^-1(x-mu), but mu is 0 in our case
  r <- solve(C)%*%X
  
  #the product of the diagonal elements Cii is the determinant of Sigma (n flops)
  log_pdf <- -0.5*n*log(2*pi) - 0.5*sum(log(diag(C))) - 0.5*t(r)%*%r
  log_pdf
}


# mvn_direct <- function(X,Sigma){
#   log_pdf <- -0.5*n*log(2*pi) - 0.5*log(det(Sigma)) - 0.5*t(X)%*%solve(Sigma)%*%(X)
#   log_pdf
# }


#algorithm that puts them together
full_sim <- function(n,mu,Sigma){
  full_result <- sim_alg1(n, mu, Sigma)
  X <- full_result[[1]]
  C <- full_result[[2]]
  pdf <- rep(NA, nrow(X))
  for(i in 1:nrow(X)){
    pdf[i] <- mvn_pdf(X[i,],C)
    print(i)
  }
  return(pdf)
}

###
### Specify values and run algorithm
###

M <- 100
n <- 1000
phi <- 0.2
mu <- rep(0,n)
Sigma <- matrix(NA, nrow=n, ncol=n)
solve(Sigma)
#creating covariance matrix
for(i in 1:nrow(Sigma)){
  for(j in 1:ncol(Sigma)){
    Sigma[i,j] <- exp(-abs(i-j)/(n*phi))
  }
}


tic()
pdf_eval <- full_sim(n,mu,Sigma)
toc()

#n=1000, just simulation
#sim_alg1, mvn_pdf .35, .34
#sim_alg2, mvn_pdf 3.22, 3.15
#sim_alg3, mvn_pdf 4.75, 4.7

##n=2000, just simulation
#sim_alg1, mvn_pdf  2.5, 2.28, 2.31
#sim_alg2, mvn_pdf  25,  24, 24
#sim_alg3, mvn_pdf  37, 38, 38

##n=3000, just simulation
#sim_alg1, mvn_pdf   7.47, 7.47, 7.6
#sim_alg2, mvn_pdf  81, 82, 82
#sim_alg3, mvn_pdf  129, 128, 131


###
### 1b) Plot the pdf values for the samples. Report the wall time for the simulation algorithm
###    and the evaluation of the pdf (jointly) using, say, system.time.
###

#plot a histogram of the values of pdf
hist(pdf_eval, main = "Histogram of 100 evaluations of mvn pdf", xlab = "Evaluation of pdf")
system.time(full_sim(n,mu,Sigma)) # 48 seconds



###
### 1c) Now that you have successfully implemented this for n = 1000, repeat the exercise 
###      for n = 2,000; n = 3,000 and so on, going up to the largest value you are able 
###      to handle. Make plots of how the computational cost scales with increasing n.
###

n <- c(1:5*1000,11:14*500)  #1:5*1000
sys_time2 <- NULL

for(k in n){
  print(paste(k, "******************************************************"))
  
  M <- 100
  phi <- 0.2
  mu <- rep(0,k)
  Sigma <- matrix(NA, nrow=k, ncol=k)
  
  #creating covariance matrix
  for(i in 1:nrow(Sigma)){
    for(j in 1:ncol(Sigma)){
      Sigma[i,j] <- exp(-abs(i-j)/(k*phi))
    }
  }
  
  print(paste("Done with making matrix"))
  
  tic()
  pdf_eval <- full_sim(k,mu,Sigma)
  time <- toc()
  time <- time$toc-time$tic
  sys_time2 <- c(sys_time2, time)
  save(sys_time2, file= "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/homework_1/data/prob_1c_sys_t2.Rdata")
}

prob_1 <- data.frame(n, sys_time2)
colnames(prob_1) <- c("n", "system_time")
prob_1$scale <- "actual_comp_time"

prob1_b <- cbind(m, m^3)
colnames(prob1_b) <- c("n", "system_time")
prob1_b$scale <- "m^3"

prob_1_full <- rbind(prob_1, prob1_b)
prob_1_full$n <- as.numeric(as.character(prob_1_full$n))
prob_1_full$system_time <- as.numeric(as.character(prob_1_full$system_time))

ggplot(prob_1_full) + geom_line(aes(x=n, y = system_time, color = scale), size = 1.2)+
  labs(title = "Computational Cost scaling with n", y = "wall time")



###
### 1d) Write out the computational complexity of your algorithm. Show the details of your 
###      calculation. If you have used an R function (e.g. eigen), you need to find out the 
###      cost of that algorithm by looking through the documentation.
###

#This is included in the report

###
### 1e) Plot the theoretical cost of your algorithm as you increase n. Do this in terms of 
###     flops (floating point operations). How does the computation scale here when compared
###     to the plot you produced above? (Try to make the plots easy to compare.) Do you have 
###     any explanation for differences between this plot and the above plot? Note: the 
###     discussion becomes more interesting as n gets large.
###

prob_1_flops <- function(n,M){
  flops_sim <- n^3/3+ (2*n-1)*M*n
  flops_eval_pdf <- M*(3*n+6)
  tot_flops <- flops_sim + flops_eval_pdf
  return(tot_flops) 
}

flop_alg <- prob_1_flops(N,100)

flop_df <- data.frame(n, flop_alg)
colnames(flop_df) <- c("n", "flops")
flop_df$flops <- as.numeric(as.character(flop_df$flops))
flop_df$n <- as.numeric(as.character(flop_df$n))

#[c(1:10, 28:37),] for zoom
ggplot(flop_df) + geom_line(aes(x=n, y = flops), size = 1.2)+
  labs(title = "n vs flops for full algorithm", ylab = "flops")



###
### Problem 2: 
###    Monte Carlo methods for random matrices: Consider the following approach for generating
###    random covariance matrices of dimension m x m: Simulate m^2 N(0; 1) random variates 
###    to obtain an mxm matrix R. Obtain the positive definite covariance matrix Sigma = RR^T. 
###    Note that this is a way to simulate a matrix with a particular Wishart distribution. 
###    Read through all the parts below before beginning to work on this problem.
###
#clear workspace
rm(list=ls())

#For an example,
m <- 100
R <- matrix(rnorm(m^2,0,1), nrow = m, ncol = m)
Sigma <- R%*%t(R)
eigen_decomp <- eigen(Sigma)
lambda <- eigen_decomp$values[1:3]
d1 <- lambda[1]-lambda[2]
d2 <- lambda[2]-lambda[3]

###
###  2a)
###

#now I will do this using Monte Carlo methods
m <- c(10000)#,1000)#,10000) skip 10,000 for now
n <- 100 #Monte Carlo sample size
d1 <- rep(NA,n)
d2 <- rep(NA,n)
df <- NULL
for(j in m){
  print(paste(j, "**********************"))
  d1 <- rep(NA,n)
  d2 <- rep(NA,n)
  for(i in 1:n){
    #if(i %% 100 == 0){ print (i)}
    print(i)
    R <- matrix(rnorm(j^2,0,1), nrow = j, ncol = j)
    Sigma <- R%*%t(R)
    #eigen_decomp <- eigen(Sigma, only.values = T)
    eigen_decomp <- eigs(Sigma, 3, opts=list(retvec=FALSE))
    #smallest <- eigs(Sigma, 1, which = "LM", sigma = 0)$values
    # https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html
    lambda <- eigen_decomp$values[1:3]
    d1_new <- lambda[1]-lambda[2]
    d2_new <- lambda[2]- lambda[3]
    d1[i] <- d1_new
    d2[i] <- d2_new
  }
  df <- cbind(df, d1,d2)
}
dfb <- df

dfb2 <- as.data.frame(df)
dfb2$d1 <- as.numeric(dfb2$d1)
dfb2$d2 <- as.numeric(dfb2$d2)

#expected value of d1 and d2
mean(dfb2[,1]) 
mean(dfb2[,2]) 

#Monte Carlo Standard Errors
sd(dfb2[,1])/sqrt(n) 
sd(dfb2[,2])/sqrt(n) 

#plot approximate density plots for d1, d2
plot(density(dfb2[,1]), "Density of d1 at m = 10,000")
plot(density(dfb2[,2]), "Density of d2 at m = 10,000")

###
###  2b)
###

#What is the expected value of the smallest eigenvalue?
m <- c(100,1000)
n <- 1000 #Monte Carlo sample size
sm <- rep(NA,n)
df <- NULL
for(j in m){
  print(paste(j, "**********************"))
  sm <- rep(NA,n)
  for(i in 1:n){
    #if(i %% 100 == 0){ print (i)}
    print(i)
    R <- matrix(rnorm(j^2,0,1), nrow = j, ncol = j)
    Sigma <- R%*%t(R)
    #eigen_decomp <- eigen(Sigma, only.values = T)
    #eigen_decomp <- eigs(Sigma, 3, opts=list(retvec=FALSE))
    smallest <- eigs(Sigma, 1, which = "LM", sigma = 0,opts=list(retvec=FALSE))$values
    # https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html
    sm[i] <- smallest[1]
  }
  df <- cbind(df,sm)
}

#expected value of d1 and d2
mean(df[,1]) #m=100
mean(df[,2]) #m=1000

#Monte Carlo Standard Errors
sd(df[,1])/sqrt(n) #m=100
sd(df[,2])/sqrt(n) #m=1000

save(df, file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/homework_1/data/prob_2b.Rdata")

#plot approximate density plots for smallest eigen value
plot(density(df[,1]), "Density of smallest eigen at m = 100")
plot(density(df[,2]), "Density of smallest eigen at m = 1,000")

###
###  2c)
###

#Now study how the expected smallest eigenvalue changes as a function of m. 
#(i) Plot approximate E(lambda_m), expected value of smallest eigenvalue, versus m. 
m <- c(1:100 * 10, 2:6*500) #(up to 3,000)
exp_val_vec <- NULL
n <- 50 #sample size
new_dat <- NULL
new_m <- NULL 
for(j in m){
  print(paste(j, "**********************"))
  sm <- rep(NA,n)
  for(i in 1:n){
    #if(i %% 100 == 0){ print (i)}
    print(i)
    R <- matrix(rnorm(j^2,0,1), nrow = j, ncol = j)
    Sigma <- R%*%t(R)
    #eigen_decomp <- eigen(Sigma, only.values = T)
    #eigen_decomp <- eigs(Sigma, 3, opts=list(retvec=FALSE))
    smallest <- eigs(Sigma, 1, which = "LM", sigma = 0,opts=list(retvec=FALSE))$values
    # https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html
    sm[i] <- smallest[1]
  }
  exp_val <- mean(sm)
  exp_val_vec <- c(exp_val_vec, exp_val)
}

prob_2 <- data.frame(m, exp_val_vec)
colnames(prob_2) <- c("m", "expected_value")
prob_2$m <- as.numeric(as.character(prob_2$m))
prob_2$expected_value <- as.numeric(as.character(prob_2$expected_value))


ggplot(prob_2) + geom_line(aes(x=m, y = expected_value), size = 1.2)+
  labs(title = "m vs expected value of smallest eigenvector", ylab = "expected value")


#(ii) Report your grid of m values. (iii) What is the largest m value for which you 
#     approximated the expectation?


###
### Problem 3
###

# Compute the inverse of the matrix using two different algorithms: 
#    (i) directly by using the solve function 
#    (ii) using the Sherman-Morrison-Woodbury identity discussed in class.


# (a) Plot the CPU time versus N for algorithm 1 and algorithm 2 
#     (you will have to determine the grid and how large you can go
#     with N)

#clear workspace
rm(list=ls())

#defining variables and data structures
sigma <- 0.2
M <- 10
N <- c((1:10)*10,(2:10)*100, (3:10)*500)
a1_time <- NULL
a2_time <- NULL

#algorithm 1, based on flops
alg_1 <- function(mat){
  inv_1 <- solve(mat)     #N^3/3 flops
  inv_1
}

#algorithm 2, based on SMW
alg_2 <- function(A, U, C, V){                    #FLOPS
  A_inv <- solve(A)                               # N
  part_1 <- A_inv%*%U                             # NxM  
  part_2 <- solve((solve(C) + V%*%A_inv%*%U))     # M, MxN, M^2, (2N-1)M^2, M^3/3  
  part_3 <- V%*%A_inv                             # MxN  
  inv_2 <- A_inv - part_1%*%part_2%*%part_3       # N^2, (2M-1)NM, (2M-1)NM      
  inv_2
  
  #inv_2 <- solve(A) - solve(A)%*%U%*%solve((solve(C) + V%*%solve(A)%*%U))%*%V%*%solve(A)
}


#using both algorithms for many values of N
for(i in N){
  print(i)
  #i <- N[1]
  norm_vec <- rnorm(i*M, 0,1)
  K <- matrix(norm_vec, nrow = i, ncol = M)
  
  A <- diag(sigma, i)
  Sig <- A + K%*%t(K)
  
  C <- diag(1, M)
  U <- K
  V <- t(K)
  
  #sum(A+ U%*%C%*%V == Sig) == nrow(Sig)*ncol(Sig) #these are the same matrix
  
  alg1_time <- system.time(alg_1(Sig))[3]
  a1_time <- c(a1_time, alg1_time)
  
  alg2_time <- system.time(alg_2(A, U, C, V))[3]
  a2_time <- c(a2_time, alg2_time)
}

alg_df <- cbind(c(N,N),c(rep("alg_1 (solve)", length(N)),rep("alg_2 (SMW)", length(N))), c(a1_time, a2_time))
colnames(alg_df) <- c("N", "algorithm", "time")
alg_df <- as.data.frame(alg_df)
alg_df$time <- as.numeric(as.character(alg_df$time))
alg_df$N <- as.numeric(as.character(alg_df$N))

#for zoom: [c(1:15, 30:43),]
ggplot(alg_df) + geom_line(aes(x=N, y = time, col = algorithm), size = 1.2)+
  labs(title = "N vs CPU for algorithm 1 and 2, zoomed in", ylab = "CPU time")

# (b) a plot of floating point operations (flops) versus N for algorithm 1 
#     and algorithm 2

#If the solve function uses cholesky decomposition, then it has n^3/3 flops
#    So, for algorithm 1, we have N^3/3 flops
#    
#    For algorithm 2, it is more complicated. I will break it up into parts.
#         solve(A) -- N^2 because diagonal
#         solve(C) -- M^2 because diagonal
#         solve(A) - solve(A)%*%U%*%solve((solve(C) + V%*%solve(A)%*%U))%*%V%*%solve(A)
#          N      N^2   N  N^2 N^2   M^2    M      N^2     N   N^2   M^2  M^2   N
#          =9N^2 + 4M^2

flops_alg_1 <- function(N){
  flops <- N^3/3
  flops
}

flops_alg_2 <- function(N, M){
  flops <- N^2 + M^2 + 3*N*M + M^2 + N^2 + M^3/3 + (2*N-1)*M^2 + 6*(2*M-1)*N^2 + (2*M-1)*N*M 
  #  N + M + 3*N*M + M^2 + N^2 + M^3/3 + (2*N-1)*M^2 + (2*M-1)*N^2 + (2*M-1)*N*M 
  flops
}


flop_alg_1 <- flops_alg_1(N)
flop_alg_2 <- flops_alg_2(N,M)

flop_df <- cbind(c(N,N),c(rep("alg_1 (solve)", length(N)),rep("alg_2 (SMW)", length(N))), c(flop_alg_1, flop_alg_2))
colnames(flop_df) <- c("N", "algorithm", "flops")
flop_df <- as.data.frame(flop_df)
flop_df$flops <- as.numeric(as.character(flop_df$flops))
flop_df$N <- as.numeric(as.character(flop_df$N))

#[c(1:10, 28:37),] for zoom
ggplot(flop_df) + geom_line(aes(x=N, y = flops, col = algorithm), size = 1.2)+
  labs(title = "N vs flops for algorithm 1 and 2", ylab = "flops")



###
### Problem 4
###

#clear workspace
rm(list =ls())

#create variables/starting values

prob_4 <- matrix(NA, nrow=4, ncol =3)
length_ci <- matrix(NA, nrow=4, ncol =3)
lambda <- c(2,2,50,50)
n <- c(10,200,10,200)
n_mc <- 10000000
sum_1 <- NULL
sum_2 <- NULL
sum_3 <- NULL

#build table for report, using n_mc MC iterations
for(i in 1:4){
  sum_1 <- NULL; ci1_f <- NULL; len_ci1 <- NULL
  sum_2 <- NULL; ci2_f <- NULL; len_ci2 <- NULL
  sum_3 <- NULL; ci3_f <- NULL; len_ci3 <- NULL
  print(i)
  
  #i = 3
  for(j in 1:n_mc){
    if(j %% 100000 == 0){ print (j)}
    
    #simulating poisson
    samp <- rpois(n= n[i], lambda = lambda[i])
    
    #taking sample mean and sample standard deviation
    est_mean <- mean(samp)
    est_sd <- sd(samp)
    
    #creating the 3 confidence intervals
    ci_1 <- cbind(est_mean - 1.96*sqrt(est_mean)/sqrt(n[i]), est_mean + 1.96*est_mean/sqrt(n[i]))
    ci_2 <- cbind(est_mean - 1.96*est_sd/sqrt(n[i]), est_mean + 1.96*est_sd/sqrt(n[i]))
    ci_3 <- cbind(qpois(p=0.025, lambda = est_mean), qpois(p=0.975, lambda = est_mean))
    
    #indicator if the true lambda is inside the confidence interval
    ind1 <- c(lambda[i] > ci_1[1] & lambda[i] < ci_1[2])
    sum_1 <- sum(sum_1,ind1)
    
    #length of confidence interval
    len_ci1[j] <- ci_1[2] - ci_1[1]
    
    #if the true lambda was outside of the confidence interval, save this confidence interval
    if(ind1 == FALSE){
      ci1_f <- rbind(ci1_f, ci_1)
    }
    
    
    #repeating for the second and third confidence interval
    ind2 <- c(lambda[i] > ci_2[1] & lambda[i] < ci_2[2])
    sum_2 <- sum(sum_2,ind2)
    len_ci2[j] <- ci_2[2] - ci_2[1]
    
    if(ind2 == FALSE){
      ci2_f <- rbind(ci2_f, ci_2)
    }
    
    ind3 <- c(lambda[i] > ci_3[1] & lambda[i] < ci_3[2])
    sum_3 <- sum(sum_3,ind3)
    len_ci3[j] <- ci_3[2] - ci_3[1]
    
    if(ind3 == FALSE){
      ci3_f <- rbind(ci3_f, ci_3)
    }
    
  }
  #can do histograms to see distribution
  #hist(len_ci1, main = expression(paste("CI with sample mean, " lambda " = 2")))
  #hist(len_ci2, main = expression(paste("CI with sample std dev, " lambda " = 2")))
  #hist(len_ci1, main = expression(paste("CI with sample mean, " lambda " = 50")))
  #hist(len_ci2, main = expression(paste("CI with sample std dev, " lambda " = 50")))
  
  #storing coverate (p)
  cov_1 <- sum_1/n_mc
  cov_2 <- sum_2/n_mc
  cov_3 <- sum_3/n_mc
  
  #storing standard error sqrt((p(1-p)/n))
  se_1 <- sqrt((cov_1*(1-cov_1))/n_mc)
  se_2 <- sqrt((cov_2*(1-cov_2))/n_mc)
  se_3 <- sqrt((cov_3*(1-cov_3))/n_mc)
  
  #saving the row of the data frame on coverage and standard error
  row_i <- c(paste(cov_1, round(se_1,4), sep = ","), paste(cov_2, round(se_2,4), sep = ","), paste(cov_3, round(se_3,4), sep = ","))
  prob_4[i,] <- row_i
  
  #saving row of the dataframe on length of confidence intervals
  len_i <- c(mean(len_ci1), mean(len_ci2), mean(len_ci3))
  length_ci[i,] <- len_i
}

xtable(prob_4)


###
### Problem 5
###
# Let X1,...,Xn iid from Gamma(alpha; beta) distribution. Consider Bayesian
# inference for alpha; beta with prior distributions a~ N(0; 3), beta~ N(0; 3). Use a
# Laplace approximation (as discussed in class) to approximate the posterior
# expectation of alpha, beta.

#clear workspace
rm(list=ls())

#load data
prob_5_dat <- read.table("C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/homework_1/data/prob_5_dat.txt", quote="\"", comment.char="")
y <- c((t(as.matrix(prob_5_dat))))
n <- length(y)

set.seed(123)

# Find unnormalized log posterior of the model
#    given parameter vector pparam and vector of data points y


tic()
denom <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  denom <- log_post
}



num <- function(param, y){
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  num <- log_post + log(param["beta"])# + log(param["beta"])
}


#give a set of initial values
initial_values <- c(alpha = 3, beta = 0.5)

#fnscale is -1 so that it maximizes the log posterior likelihood
opt_fit_den <- optim(initial_values, denom, control = list(fnscale = -1), y = y, hessian = TRUE)
opt_fit_num <- optim(initial_values, num, control = list(fnscale = -1), y = y, hessian = TRUE)

full_den <- exp(opt_fit_den$value)*(det(solve(opt_fit_den$hessian))^(-1/2))
full_num <- exp(opt_fit_num$value)*(det(solve(opt_fit_num$hessian))^(-1/2))

exp_est <- full_num/full_den
toc()