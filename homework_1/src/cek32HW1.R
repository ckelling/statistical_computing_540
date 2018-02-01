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
M <- 100
n <- 1000
phi <- 0.9
mu <- rep(0,n)
Sigma <- matrix(NA, nrow=n, ncol=n)

#creating covariance matrix
for(i in 1:nrow(Sigma)){
  for(j in 1:ncol(Sigma)){
    Sigma[i,j] <- 2*exp(-abs(i-j)/phi)
  }
}

#simulating 100 draws, through R function
sim_mvnorm <- mvrnorm(M, mu, Sigma)

#by hand
X <- matrix(rnorm(n*M), n)
L <- chol(Sigma)
sum(L%*%t(L) == Sigma) #not exact, but close
#create multivariate normal
#https://stats.stackexchange.com/questions/12953/generating-values-from-a-multivariate-gaussian-distribution
#https://en.wikipedia.org/wiki/Cholesky_decomposition
mv_n <- mu + L%*%X

sim_alg <- function(M){
  X <- matrix(rnorm(n*M), n)
  L <- chol(Sigma)
  mv_n <- mu + L%*%X
  return(mv_n)
}

# e_decomp <- eigen(Sigma)
# e_values <- e_decomp$values
# X <- matrix(rnorm(n * M), M)
# X <- mu + e_decomp$vectors %*% diag(sqrt(e_values), n) %*% t(X)

#Write your own code to also evaluate the multivariate normal pdf at each simulated value.
#   I will evaluate the log pdf and then convert it back
#   https://www.cs.cmu.edu/~epxing/Class/10701-08s/recitation/gaussian.pdf
#   mu is 0 in this case so I will take it out

#Calculate pdf for each X
log_pdf_vec <- rep(NA, M)
for(i in 1:M){
  X1 <- mv_n[,i]
  #X1 <- sim_mvnorm[1,]
  log_pdf <- -0.5*n*log(2*pi) - 0.5*log(det(Sigma)) - 0.5*t(X1)%*%solve(Sigma)%*%(X1)
  log_pdf_vec[i] <- exp(log_pdf)
}


###
### 1b) Plot the pdf values for the samples. Report the wall time for the simulation algorithm
###    and the evaluation of the pdf (jointly) using, say, system.time.
###

#plot a histogram of the values of pdf
hist(log_pdf_vec)

#simulation algorithm
full_sim <- function(M){
  mv_n <- sim_alg(M)
  log_pdf_vec <- rep(NA, M)
  for(i in 1:M){
    X1 <- mv_n[,i]
    #X1 <- sim_mvnorm[1,]
    log_pdf <- -0.5*n*log(2*pi) - 0.5*log(det(Sigma)) - 0.5*t(X1)%*%solve(Sigma)%*%(X1)
    log_pdf_vec[i] <- exp(log_pdf)
  }
  return(log_pdf_vec)
}

system.time(full_sim(100)))

###
### 1c) Now that you have successfully implemented this for n = 1000, repeat the exercise 
###      for n = 2,000; n = 3,000 and so on, going up to the largest value you are able 
###      to handle. Make plots of how the computational cost scales with increasing n.
###

M <- 1:3*1000
sys_time <- NULL
for(i in M){
  sys_t <- system.time(full_sim(M)))
  sys_time <- c(sys_time, sys_t)
}


###
### 1d) Write out the computational complexity of your algorithm. Show the details of your 
###      calculation. If you have used an R function (e.g. eigen), you need to find out the 
###      cost of that algorithm by looking through the documentation.
###

###
### 1e) Plot the theoretical cost of your algorithm as you increase n. Do this in terms of 
###     flops (floating point operations). How does the computation scale here when compared
###     to the plot you produced above? (Try to make the plots easy to compare.) Do you have 
###     any explanation for differences between this plot and the above plot? Note: the 
###     discussion becomes more interesting as n gets large.
###


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
m <- c(1000)#,1000)#,10000) skip 10,000 for now
n <- 1000 #Monte Carlo sample size
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
    d1_new <- lambda[3]-lambda[2]
    d2_new <- lambda[2]-lambda[1]
    d1[i] <- d1_new
    d2[i] <- d2_new
  }
  df <- cbind(df, d1,d2)
}
dfb <- df
#expected value of d1 and d2
mean(d1) # m=100: 9.69232    , m=1,000: 9.591571    , m=10,000: 9.537369
mean(d2) # m=100: 6.271687   , m=1,000: 6.284819    , m=10,000: 6.324104

#Monte Carlo Standard Errors
sd(d1)/sqrt(n)  # m=100: 5.546734   , m=1,000: 5.392687    , m=10,000: 5.470581
sd(d2)/sqrt(n)  # m=100: 3.493466   , m=1,000: 3.499043    , m=10,000: 3.52086

#plot approximate density plots for d1, d2
plot(density(d1), "Density of d1 at m = 100")
plot(density(d2), "Density of d2 at m = 100")

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
m <- c(1:100 * 10,3:6*500) #1:100 * 10,
test <- 1:100 * 10
#exp_val_vec2 <- exp_val_vec
exp_val_vec <- NULL
n <- 100 #sample size

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

###

###

prob_2 <- data.frame(m, exp_val_vec)
colnames(prob_2) <- c("m", "expected_value")
prob_2$m <- as.numeric(as.character(prob_2$m))
prob_2$expected_value <- as.numeric(as.character(prob_2$expected_value))

ggplot(prob_2) + geom_line(aes(x=m, y = expected_value), size = 1.2)+
  labs(title = "m vs expected value of smallest eigenvector", ylab = "expected value")

#(ii) Report your grid of m values. (iii) What is the largest m value for which you 
#     approximated the expectation?

exp_val_vec_2 <- c(0.0733591536, 0.0379156155, 0.0198245616, 0.0180753816, 0.0103755242, 0.0089777114, 0.0096156783, 0.0102288524, 0.0069491758, 0.0062841157, 0.0053903295, 0.0063654489, 0.0050404595, 0.0052720499, 0.0054420552, 0.0035475238, 0.0039151589, 0.0041045611, 0.0040739881, 0.0037872857, 0.0037378989, 0.0035142535, 0.0026135692, 0.0027201996, 0.0027336333, 0.0027757366, 0.0022619196, 0.0019878161, 0.0016583414, 0.0021004911, 0.0017371396, 0.0017891204, 0.0022129804, 0.0017180198, 0.0024197094, 0.0022395280, 0.0015340470, 0.0012885595, 0.0018150200, 0.0018829293, 0.0017238653, 0.0021366696, 0.0019007053, 0.0014812419, 0.0016249169, 0.0014746681, 0.0015488306, 0.0015422615, 0.0014832698, 0.0013745350, 0.0010942752, 0.0012683243, 0.0012176643, 0.0010465356, 0.0013874523, 0.0009151464, 0.0009779324, 0.0011047335, 0.0011135144, 0.0011581533 0.0013261632 0.0011842986 0.0009418504 0.0008637009 0.0012042640 0.0014792420 0.0009825543 0.0009585732 0.0008142311 0.0008808660 0.0009690631 0.0012030432 0.0008774910 0.0009821160 0.0007181296 0.0009557500 0.0009308187 0.0011616895 0.0009327619 0.0007620769 0.0010152202 0.0009325263 0.0007143623 0.0006568169 0.0007831074 0.0008342921 0.0007174419 0.0009398580 0.0008389230 0.0007639229 0.0006491160 0.0007129268 0.0006939538 0.0007379282 0.0007542634 0.0006987829 0.0008465385 0.0006206839 0.0006725681 0.0008509256


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
alg_2 <- function(A, U, C, V){                    #FLOPS, N 5, M 10
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
  flops <- N + 3*N*M + M + M^2 + (2*N-1)*M^2 + M^3/3 + N^2 + 2*(2*M-1)*N*M    
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
n_mc <- 1000000
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
    if(j %% 10000 == 0){ print (j)}
    
    #simulating poisson
    samp <- rpois(n= n[i], lambda = lambda[i])
  
    #taking sample mean and sample standard deviation
    est_mean <- mean(samp)
    est_sd <- sd(samp)
  
    #creating the 3 confidence intervals
    ci_1 <- cbind(est_mean - 1.96*est_mean/sqrt(n[i]), est_mean + 1.96*est_mean/sqrt(n[i]))
    ci_2 <- cbind(est_mean - 1.96*est_sd/sqrt(n[i]), est_mean + 1.96*est_sd/sqrt(n[i]))
    ci_3 <- cbind(qpois(p=0.025, lambda = lambda[i]), qpois(p=0.975, lambda = lambda[i]))
    
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

set.seed(123)

# Find unnormalized log posterior of the model
#    given parameter vector pparam and vector of data points y

log_posterior <- function(param, y) {
  log_lik <- sum(dgamma(y, param["alpha"], param["beta"], log = T))  # the log likelihood
  log_post <- log_lik + dnorm(param["alpha"], 0, 3, log = T) + dnorm(param["beta"], 0, 3, log = T)
  log_post
}

#give a set of initial values
initial_values <- c(alpha = 4, beta = 4)

#fnscale is -1 so that it maximizes the log posterior likelihood
opt_fit <- optim(initial_values, log_posterior, control = list(fnscale = -1), y = y)

#posterior estimates
post_est <- opt_fit$par
post_est

