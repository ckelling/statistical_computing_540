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

###
### Problem 1
###

# Write computer code to simulate M = 100 draws from a multivariate normal
#     distribution with mean 0 (vector of 0s) and covariance matrix
#     as above for n = 1000.
#     Sigma_{ij} = exp(-|i-j|/phi),
M <- 100
n <- 1000
phi <- 0.2
mu <- rep(0,n)
Sigma <- matrix(NA, nrow=n, ncol=n)

#creating covariance matrix
for(i in 1:nrow(Sigma)){
  for(j in 1:ncol(Sigma)){
    Sigma[i,j] <- exp(-abs(i-j)/phi)
  }
}

#simulating 100 draws
sim_mvnorm <- mvrnorm(M, mu, Sigma)

#by hand
X <- matrix(rnorm(n*M), n)
L <- chol(Sigma)
sum(L%*%t(L) == Sigma) #not exact, but close
#create multivariate normal
#https://stats.stackexchange.com/questions/12953/generating-values-from-a-multivariate-gaussian-distribution
#https://en.wikipedia.org/wiki/Cholesky_decomposition
mv_n <- mu + L%*%X

# e_decomp <- eigen(Sigma)
# e_values <- e_decomp$values
# X <- matrix(rnorm(n * M), M)
# X <- mu + e_decomp$vectors %*% diag(sqrt(e_values), n) %*% t(X)

#Write your own code to also evaluate the multivariate normal pdf at each simulated value.
#   I will evaluate the log pdf and then convert it back
#   https://www.cs.cmu.edu/~epxing/Class/10701-08s/recitation/gaussian.pdf
#   mu is 0 in this case so I will take it out
log_pdf <- -0.5*n*log(2*pi) - 0.5*log(det(Sigma)) - 0.5*t(mv_n)%*%solve(Sigma)%*%(mv_n)
pdf <- exp(log_pdf)

mvrnorm
###
### 1b) Plot the pdf values for the samples. Report the wall time for the simulation algorithm
###    and the evaluation of the pdf (jointly) using, say, system.time.
###

###
### 1c) Now that you have successfully implemented this for n = 1000, repeat the exercise 
###      for n = 2,000; n = 3,000 and so on, going up to the largest value you are able 
###      to handle. Make plots of how the computational cost scales with increasing n.
###

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
tic()
m <- 10000
n <- 10000 #sample size
d1 <- NULL
d2 <- NULL
for(i in 1:n){
  m <- 10
  R <- matrix(rnorm(m^2,0,1), nrow = m, ncol = m)
  Sigma <- R%*%t(R)
  eigen_decomp <- eigen(Sigma)
  lambda <- eigen_decomp$values[(length(eigen_decomp$values)-2):length(eigen_decomp$values)]
  d1_new <- lambda[3]-lambda[2]
  d2_new <- lambda[2]-lambda[1]
  d1 <- rbind(d1,d1_new)
  d2 <- rbind(d2,d2_new)
}
toc()

#expected value of d1 and d2
mean(d1) # m=100: 9.69232    , m=1,000: 9.591571    , m=10,000: 9.537369
mean(d2) # m=100: 6.271687   , m=1,000: 6.284819    , m=10,000: 6.324104

#Monte Carlo Standard Errors
sd(d1)  # m=100: 5.546734   , m=1,000: 5.392687    , m=10,000: 5.470581
sd(d2)  # m=100: 3.493466   , m=1,000: 3.499043    , m=10,000: 3.52086

#plot approximate density plots for d1, d2
plot(density(d1), "Density of d1 at m = 100")
plot(density(d2), "Density of d2 at m = 100")

###
###  2b)
###

#What is the expected value of the smallest eigenvalue?
tic()
m <- 10000
n <- 10000 #sample size
sm_eigen <- NULL
for(i in 1:n){
  R <- matrix(rnorm(m^2,0,1), nrow = m, ncol = m)
  Sigma <- R%*%t(R)
  eigen_decomp <- eigen(Sigma)
  sm_new <- eigen_decomp$values[length(eigen_decomp$values)]
  sm_eigen <- c(sm_eigen, sm_new)
}
toc()

#expected value of d1 and d2
mean(sm_eigen) # m=100: 9.69232    , m=1,000: 9.591571    , m=10,000: 9.537369

#Monte Carlo Standard Errors
sd(sm_eigen)  # m=100: 5.546734   , m=1,000: 5.392687    , m=10,000: 5.470581

#plot approximate density plots for d1, d2
plot(density(sm_eigen), "Density of smallest eigen at m = 100")

###
###  2c)
###

#Now study how the expected smallest eigenvalue changes as a function of m. 
#(i) Plot approximate E(lambda1), expected value of smallest eigenvalue, versus m. 
a <- 1:100
b <- 1
m <- a * b
exp_val_vec <- NULL
tic()
for(j in m){
  n <- 10000 #sample size
  sm_eigen <- NULL
  for(i in 1:n){
    R <- matrix(rnorm(j^2,0,1), nrow = j, ncol = j)
    Sigma <- R%*%t(R)
    eigen_decomp <- eigen(Sigma)
    sm_new <- eigen_decomp$values[length(eigen_decomp$values)]
    sm_eigen <- c(sm_eigen, sm_new)
  }
  exp_val <- mean(sm_eigen)
  exp_val_vec <- c(exp_val_vec, exp_val)
  print(m[j])
}
toc()
plot(x= m,y= exp_val_vec)

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

# (b) a plot of floating point operations (flops) versus N for algorithm 1 
#     and algorithm 2

#     LOOK AT ASSIGNMENT FOR C


###
### Problem 4
###

prob_4 <- matrix(NA, nrow=4, ncol =6)
lambda <- c(2,2,50,50)
n <- c(10,200,10,200)

for(i in 1:4){
  #i=1
  samp <- rpois(n= n[i], lambda = lambda[i])
  
  est_mean <- mean(samp)
  est_sd <- sd(samp)
  
  ci_1 <- cbind(est_mean - 1.96*est_mean/sqrt(n[i]), est_mean + 1.96*est_mean/sqrt(n[i]))
  ci_2 <- cbind(est_mean - 1.96*est_sd/sqrt(n[i]), est_mean + 1.96*est_sd/sqrt(n[i]))
  ci_3 <- cbind(qpois(p=0.025, lambda = lambda[i]), qpois(p=0.975, lambda = lambda[i]))
  
  total_ci <- cbind(ci_1, ci_2, ci_3)
  prob_4[i,] <- total_ci
}

###
### Problem 5
###
# Let X1; : : : ;Xn iid from Gamma(alpha; beta) distribution. Consider Bayesian
# inference for alpha; beta with prior distributions a~ N(0; 3), beta~ N(0; 3). Use a
# Laplace approximation (as discussed in class) to approximate the posterior
# expectation of alpha, beta.
