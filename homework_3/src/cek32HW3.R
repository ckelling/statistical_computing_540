###
### Claire Kelling
### STAT 540
### Homework #3
###
### Created 3/27/2018 for the third assigment, due 4/5/2018
### The exercises mainly focus on optimization algorithms
### 

library(data.table)

###
### Problem 1
###

# Find the MLE for a logistic regression model.
# Write a Newton-Raphon algorithm to find the MLE (\hat{\beta_1}; \hat{\beta_2}). 
# The data for this problem are available here: http://personal.psu.edu/muh10/540/data/logReg.dat

#clear workspace
rm(list=ls())

#load data
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")

##########################################
## Newton-Raphson
##########################################
#Find the MLE for a logistic regression model.} For $i = 1,...,n$
#  $$ Y_i \sim  Ber(p_i(\theta)); where  \theta= (\beta_1, \beta_2)$$
#    $$p_i(\theta) = exp(\beta_1X_{1i} + \beta_2X_{2i})/(1 + exp(\beta_1X_{1i} + \beta_2X_{2i}))$$ 
# where each $Y_i$ is a binary response corresponding to predictors $X_{1i},X_{2i}$. 
# Write a Newton-Raphson algorithm to find the MLE $( \hat{\beta}_1, \hat{\beta}_2)$. 

set.seed(1)
N=nrow(y)
x1=prob1$x1
x2=prob1$x2
y=prob1$y

plot(prob1$x1, prob1$y, main="x1 vs y")
plot(prob1$x2, prob1$y, main="x2 vs y")

#initial values
beta0=3
beta1=-1
beta2=1
sigmasq=1

f.lr.p <- function(X, beta){
  # compute vector p of probabilities for logistic regression with logit link
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

f.lr.l <- function(y, m, p){
  # bernoulli log likelihood function
  # input:   vectors: y = counts; p = probabilities
  # output: log-likelihood l, a scalar
  one_vec <- rep(1,nrow(prob1))
  l <- t(y) %*% log(p) + t(one_vec - y) %*% log(1 - p)
  return(l)
}
  
f.lr.FS <- function(X, y, m, beta.1, eps1 = 1e-6, eps2 = 1e-7, maxit = 50){
  # Fisher's scoring:
  #   X        = n-by-(r+1) design matrix
  #   y        = n-by-1 vector of success counts
  #   m        = n-by-1 vector of sample sizes
  #   beta.1   = (r+1)-by-1 vector of starting values for regression est
  # Iteration controlled by:
  #   eps1  = absolute convergence criterion for beta
  #   eps2  = absolute convergence criterion for log-likelihood
  #   maxit = maximum allowable number of iterations
  # Output:
  #   out   = list containing:
  #     beta.MLE  = beta MLE
  #     NR.hist   = iteration history of convergence differences
  #     beta.hist = iteration history of beta
  #     beta.cov  = beta covariance matrix (inverse Fisher's information matrix at MLE)
  #     note      = convergence note
  beta.2 <- rep(-Inf, length(beta.1)) # init beta.2
  
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  llike.1 <- f.lr.l(y, m, f.lr.p(X, beta.1)) # update loglikelihood
  llike.2 <- f.lr.l(y, m, f.lr.p(X, beta.2)) # update loglikelihood

  
  diff.like <- abs(llike.1 - llike.2) # diff
  if (is.nan(diff.like)){
    diff.like <- 1e9
  }
  i <- 1                                   # initial iteration index
  alpha.step <- seq(-1, 2, by = 0.1)[-11]  # line search step sizes, excluding 0
  
  
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1, step.size = 1) # iteration history
  beta.hist <- matrix(beta.1, nrow = 1)
  
  while((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2)){
    i <- i + 1                       # increment iteration
    # update beta
    beta.2 <- beta.1                 # old guess is current guess
    mu.2 <-  f.lr.p(X, beta.2)    # p is mean
    # variance matrix
    v.2  <- diag(as.vector(f.lr.p(X, beta.2) * (1 - f.lr.p(X, beta.2))))
    score.2 <- t(X) %*% (y - mu.2)   # score function
    # this increment version inverts the information matrix
    # Iinv.2 <- solve(t(X) %*% v.2 %*% X)  # Inverse information matrix
    # increm <- Iinv.2 %*% score.2     # increment, solve() is inverse
    # this increment version solves for (beta.2-beta.1) without inverting Information
    increm <- solve(t(X) %*% v.2 %*% X, score.2) # solve for increment
    # line search for improved step size
    llike.alpha.step <- rep(NA, length(alpha.step)) # init llike for line search
    for (i.alpha.step in 1:length(alpha.step)){
      llike.alpha.step[i.alpha.step] <- f.lr.l(y, f.lr.p(X, beta.2 + alpha.step[i.alpha.step] * increm))
    }
    
    
    # step size index for max increase in log-likelihood (if tie, [1] takes first)
    ind.max.alpha.step <- which(llike.alpha.step == max(llike.alpha.step))[1]
    beta.1 <- beta.2 + alpha.step[ind.max.alpha.step] * increm   # update beta
    diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
    llike.2 <- llike.1               # age likelihood value
    llike.1 <- f.lr.l(y, f.lr.p(X, beta.1)) # update loglikelihood
    diff.like <- abs(llike.1 - llike.2) # diff
    
    
    # iteration history
    NR.hist   <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1, alpha.step[ind.max.alpha.step]))
    beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
  }
  # prepare output
  out <- list()
  out$beta.MLE  <- beta.1
  out$iter      <- i - 1
  out$NR.hist   <- NR.hist
  out$beta.hist <- beta.hist
  v.1  <- diag(as.vector(m * f.lr.p(X, beta.1) * (1 - f.lr.p(X, beta.1))))
  Iinv.1 <- solve(t(X) %*% v.1 %*% X)  # Inverse information matrix
  out$beta.cov  <- Iinv.1
  
  if(!(diff.beta > eps1) & !(diff.like > eps2)){
    out$note <- paste("Absolute convergence of", eps1, "for betas and",
                      eps2, "for log-likelihood satisfied")
  }
  if(i > maxit){
    out$note <- paste("Exceeded max iterations of ", maxit)
  }
  return(out)
}

# (a) Provide details about how you set up the algorithm, including any
#  analytical quantities you derived. You should also include the stopping criteria, 
#  how you obtained initial values, how you tuned sn (step size/learning rate) and other 
#  relevant information, e.g. did you run the algorithm using multiple initial values?

n <- nrow(prob1)
y <- prob1$y
x1 <- prob1$x1
x2 <- prob1$x2
# quadratic model
X <- matrix(c(x1, x2), nrow = n)
colnames(X) <- c("x1", "x2")
r <- ncol(X)

# initial beta vector
beta.1 <- rep(0, r)
# fit betas using our Fisher Scoring function
out <- f.lr.FS(X, y, beta.1)
out



set.seed(101)
n.obs <- 115
n.zero <- 60
n.pred <- 30
y <- c(rep(0,n.zero),rep(1,n.obs-n.zero))
X <- sample(c(0,1),size=n.pred*n.obs,replace=TRUE)
Z <- t(matrix(X,ncol=n.obs))

g1 <- glm(y~.-1,data.frame(y,Z),family="binomial")  

NRfit <- function(y,X,start,n.iter=100,tol=1e-4,verbose=TRUE) {
  ## used X rather than Z just because it's more standard notation
  n.pred <- ncol(X)
  B <-  matrix(NA,ncol=n.iter,
               nrow=n.pred)
  B[,1] <- start
  for (i in 2:n.iter) {
    if (verbose) cat(i,"\n")
    p <- plogis(X %*% B[,i-1])
    v.2 <- diag(c(p*(1-p)))
    score.2 <- t(X) %*% (y - p) # score function
    increm <- solve(t(X) %*% v.2 %*% X) 
    B[,i] <- B[,i-1]+increm%*%score.2
    if (all(abs(B[,i]-B[,i-1]) < tol)) return(B)
  }
  B
}

matplot(res1 <- t(NRfit(y,Z,start=c(0,0))))
matplot(res2 <- t(NRfit(y,Z,start=rep(0,ncol(Z)))))
all.equal(res2[6,],unname(coef(g1))) ## TRUE

# (b) Provide your pseudocode for the main algorithm. You can use notation defined in class 
#  (and specified above) in order to make the pseudocode succinct and easy to read. 

# (c) Provide the MLE as well as estimate of the standard error, along
#  with asymptotic 95% confidence intervals.

# (d) Report the total computing time taken by your algorithm (e.g. you
#  could use system.time). If you have 2-3 different versions of the
#  algorithm, report timings for each version.

# (e) Provide the computational cost of the algorithm (flops) per iteration
#  as a function of n. Also provide the number of iterations it took before
#  the stopping criteria was met. Of course, the number of iterations will
#  vary by initial conditions so you should provide at least 2-3 different
#  counts depending on where you started.

# (f) Summarize in 1-2 sentences the sensitivity of the algorithm to the
#  initial value, that is, does the algorithm work well regardless of where
#  you start or does it need to be in a neighborhood of a certain value?
#    (This is not required, but because this is a simple 2-D optimization
#     problem, you could also provide the 2-D log-likelihood surface to
#     obtain insights.)



###
### Problem 2
###

# Implement a gradient descent algorithm (work with the negative log like-lihood). 
# Answer (a)-(f) for this algorithm.

contour(bvals1,bvals2,z)
points(beta1,beta2, pch=19, col="red")
NUMIT=100
est=c(-1,-2) # initial value
sRate=1 #1.2 2, 0.5
points(est[1], est[2], pch=19, col="blue")
print(est)
for (i in 1:NUMIT){
  df1=-(2/N)*sum(x1*(y-est[1]*x1-est[2]*x2))
  df2=-(2/N)*sum(x2*(y-est[1]*x1-est[2]*x2))
  ##    df1=-2*sum(x1*(y-est[1]*x1-est[2]*x2))
  ##    df2=-2*sum(x2*(y-est[1]*x1-est[2]*x2))
  
  est[1]=est[1]-sRate*df1
  est[2]=est[2]-sRate*df2
  cat(paste("i=",i,", est=(",est[1],",",est[2],") \n"))
  points(est[1], est[2], pch=19, col="blue")
  ##    points(est[1], est[2], pch=paste(i), col="blue")
  readline("next iteration")
  ##    text(est[1],est[2],label=i,col="blue)
}

###
### Problem 3
###

# Now implement a stochastic gradient descent algorithm and repeat (a)-(f).
set.seed(1)
contour(bvals1,bvals2,z)
points(beta1,beta2, pch=19, col="red")
NUMIT=500
est=c(1,1) # initial value
##est=c(-1,-2) # initial value
sRate= 3 #2 1 0.5 5
points(est[1], est[2], pch=19, col="blue")
print(est)
for (i in 1:NUMIT){
  randInd= sample(seq(1,N), 1)
  df1=-(2/N)*sum(x1[randInd]*(y[randInd]-est[1]*x1[randInd]-est[2]*x2[randInd]))
  df2=-(2/N)*sum(x2[randInd]*(y[randInd]-est[1]*x1[randInd]-est[2]*x2[randInd]))
  
  est[1]=est[1]-sRate*df1
  est[2]=est[2]-sRate*df2
  cat(paste("i=",i,", est=(",est[1],",",est[2],") \n"))
  points(est[1], est[2], pch=19, col="blue")
  ##    points(est[1], est[2], pch=paste(i), col="blue")
  readline("next iteration")
  ##    text(est[1],est[2],label=i,col="blue)
}



###
### Problem 4
###

#  Compare the three algorithms above: Newton-Raphson, gradient descent,
# and stochastic gradient descent. Provide a 1-2 sentence summary of which
# algorithm you would recommend and why. Then provide a more detailed
# comparison, for example how stable the algorithms are with respect to
# initial values, how sensitive they are to choice of s_n, a comparison of the
# total computational cost, and the computational cost per iteration.


###
### Problem 5
###

# E-M algorithm: Return to a lightbulb lifetimes problem similar to the
# problem in a previous homework. Assume lightbulbs made by a par-
#   ticular company are independent and gamma distributed with param-
#   eters alpha; beta (parameterization is such that expected value is alpha*beta). Sup-
#   pose in an experiment m bulbs are switched on at the same time, but
# are only completely observed up to time tau . Let the lifetimes of these
# bulbs be A1;A2; : : : ;Am. However, since the bulbs are only observed
# till time tau, not all these lifetimes will be observed. Now suppose that
# at time tau, the experimenter observes the number of lightbulbs, W, still
# working at time tau , and the lifetimes of all lightbulbs that stopped work-
#   ing by tau . For convenience, denote these bulbs that stopped working by
# time tau as A1;A2; : : : ;Am????W. Hence, the missing information consists of
# the lifetimes of the lightbulbs still working at time tau , Am????W+1; : : : ;Am.
# For a particular experiment, let tau be 200 days and m = 300. The
# data on the lightbulb lifetimes for the bulbs that stopped working by
# tau are here: http://personal.psu.edu/muh10/540/data/bulbsHW3.dat
# Assume that the remaining bulbs were still working at time tau . Find the
# MLE for (alpha; beta) using the E-M algorithm. Repeat parts (a)-(f) from above
# for this algorithm as well.

times <- fread("http://personal.psu.edu/muh10/540/data/bulbsHW3.dat")