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

prob_5 <- scan("http://personal.psu.edu/muh10/540/data/bulbsHW3.dat")
data <- as.vector(prob_5)

tau <- 200 #stopping time in number of days
m <- 300
obs <- length(data[,1]) #number of observed lifetimes
W <- m-obs #number of unobserved lifetimes



#### estimate using 50 as the starting value ####
a.1 <- as.numeric(rep(NA, 1000))
b.1 <- as.numeric(rep(NA, 1000))
a.1[1] <- 2 #given starting value
b.1[1] <- 2 #given starting value


# We evaluate the first estimate outside of the loop so that a while condition can be applied later
# The Q-function is based on the log-likelihood
Q.fn <- function(param){
  a <- param[1]
  b <- param[2]
  (-m*a*log(b)) - m*lgamma(a)+(a - 1)*(sum(log(data))+log(W*(tau+a.1[1]*b.1[1])))-
    ((1/b)*(sum(data) + (W*(tau+a.1[1]*b.1[1]))))
  
  #(-m*a*log(b)) - m*lgamma(a)+(a - 1)*(sum(log(data)))-
  #  ((1/b)*(sum(data) + (W*(tau+a.1[1]*b.1[1]))))
}

test <- optim(par=c(a.1[1], b.1[1]), fn=Q.fn,
                         method="L-BFGS-B",
              lower=c(0,0), upper=c(20,20), 
              control=list(fnscale=-1))$par

# We will now optimize the parameter value until it converges
# Convergence criteria: if the difference between (k-1)th and kth estimate
# is less than 0.5, we accept the kth value

eps <- 1e-4
while ((a.1[max(which(!is.na(a.1)))])- (a.1[max(which(!is.na(a.1)))-1]) > eps|
       (b.1[max(which(!is.na(b.1)))])- (b.1[max(which(!is.na(b.1)))-1]) > eps){
  
  # We define the Q-function inside the loop to make reference to kth iterate easier
  # The Q-function is based on the log-likelihood
  Q.fn <- function(param){
    a <- param[1]
    b <- param[2]
    (-m*a*log(b)) - m*lgamma(a)+(a - 1)*
      (sum(log(data))+log(W*(tau+a.1[max(which(!is.na(a.1)))]*b.1[max(which(!is.na(b.1)))])))-
      ((1/b)*(sum(data) + (W*(tau+a.1[max(which(!is.na(a.1)))]*b.1[max(which(!is.na(b.1)))]))))
    
  }
  
  # We now optimize the value of theta
  new_fit <- optim(par=as.numeric(theta.vals.1[max(which(!is.na(theta.vals.1)))]), fn=Q.fn,
                                                            method="Brent", lower=0, upper=9999,
                                                            control=list(fnscale=-1))$par
  test <- optim(par=c(a.1[1], b.1[1]), fn=Q.fn,
                method="L-BFGS-B",
                lower=c(0,0), upper=c(20,20), 
                control=list(fnscale=-1))$par
  
  a.1[max(which(!is.na(a.1)))+1] <- new_fit[1]
  b.1[max(which(!is.na(b.1)))+1] <- new_fit[2]
}

plot(na.omit(theta.vals.1), main="Estimated theta: Starting value 50", xlab="Number of iterations", ylab="Estimated theta", col="dodgerblue4", pch=19)
# First estimate of theta
theta.hat.1 <- theta.vals.1[max(which(!is.na(theta.vals.1)))]
theta.hat.1


library(truncdist)
x <- seq( 0, 3, .1 )
pdf <- dtrunc(x, spec="gamma", a=0, b=20, shape =2, scale =2 ) #pdf of truncated gamma
test <- extrunc(spec = "gamma", a = 0, b = 3, shape = 2, scale = 9)
pdf2  <- dtrunc(x, spec="gamma", a=1, shape =2, scale =2 ) #pdf of truncated gamma
hist(pdf2)
hist(x)

Q.function <- dtrunc(x, spec="gamma", a=0, b=tau, shape =2, scale =2 )*
  dtrunc(x, spec="gamma", a=tau, shape =2, scale =2 )*dbinom(x,n,p)

library(pracma)

x <- 3
dgamma3(9, shape =1, scale = 1, thres = 5)
dgamma(3, shape =1, scale = 1)
dtrunc(9, spec="gamma", shape = 1, scale = 1, a = 5)

extrunc(spec = "gamma", b = 2, shape = 1, scale = 1)

test <- dtrunc(data,spec="gamma", shape = alpha, scale = beta, b = tau, log = F)*
        dtrunc(unobs,spec="gamma",shape=alpha.o, scale = beta.o, a = tau, log = F)*
        dbinom(n=m, p = (1-(incgam(w,alpha.o)/gamma(alpha.o))))

unobs <- rtrunc(W, spec = "gamma", shape = alpha, scale = beta, a = tau)
test3 <- dtrunc(unobs, spec = "gamma", shape = alpha, scale = beta, a = tau, log =T)
sum(test3)

Q_func <- function(param,a.o, b.o){
  alpha <- param[1]
  beta <- param[2]
  like <- sum(dtrunc(data,spec="gamma", shape = alpha, scale = beta, b = tau, log = T))+
    sum(dtrunc(unobs,spec="gamma",shape=a.o, scale = b.o, a = tau, log = T))+
    dbinom(x = W, size = m, p = (1-(incgam(W,a.o)/gamma(a.o))), log = T)
  return(like)
}

#compute optimal
a.2 <- 10
b.2 <- 10
a.o <- 10
b.o <- 10
init <- c(a.2,b.2)
params <- optim(par=init, a.o =a.2, b.o =b.2, fn=Q_func,
                method="L-BFGS-B",
                lower=c(0,0), upper=c(10,10), 
                control=list(fnscale=-1))$par


unobs_2 <- rep(tau,W)
Q_paper <- function(param){
  alpha <- param[1]
  beta <- param[2]
  log_like <- sum((alpha-1)*data - data/beta - alpha*log(beta) - lgamma(alpha))+
              W*(log(1-(gammainc(alpha,tau/beta)[1]/gamma(alpha))))
  return(log_like)
}
test <- Q_paper(alpha,beta)
gammainc(2,2)[1]*3

alpha <- 1 ; beta <- 2

init <- c(alpha,beta)
params <- optim(par=init, a.o =a.2, b.o =b.2, fn=Q_func,
                method="L-BFGS-B",
                lower=c(0,0), upper=c(10,10), 
                control=list(fnscale=-1))$par



Q_combined <- function(alpha,beta,a.o,b.o){
  
  log_like <- sum(dtrunc(data,spec="gamma", shape = alpha, scale = beta, b = tau, log = T))+
    sum(dtrunc(unobs,spec="gamma",shape=a.o, scale = b.o, a = tau, log = T))+
    dbinom(x = W, size = m, p = (1-(incgam(W,a.o)/gamma(a.o))), log = T)
  
  exp_log_like <- sum((alpha-1)*data - data/beta - alpha*log(beta) - lgamma(alpha))+
    W*(log(1-(gammainc(alpha,tau/beta)[1]/gamma(alpha))))
  
  
}



delta <- 1e-6

approx <- 1/delta * (incgam(alpha + delta, tau/beta)- incgam(alpha, tau/beta))
1/delta*(gammainc(alpha + delta, tau/beta)[2]- gammainc(alpha, tau/beta)[2])
1/delta*(incgam(alpha + delta, tau/beta)- incgam(alpha, tau/beta))


alpha <- 90

test2 <- cbind(data, test)
prod3 <- prod(test2)

#survival probability, for binomial distribution
#https://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
s_p <- 1- gammainc(x,a2)/gamma(x) # a is the
gammainc(x,a2)






###
### Coding up Paper
###

C2 <- function(alpha, beta, tau){
  num <- beta*incgam(alpha+1, tau/beta)
  denom <- incgam(alpha, tau/beta)
  e_2 <- num/denom
  return(e_2)
}


C1 <- function(alpha, beta, tau){
  #approximating derivative
  delta <- 1e-6
  approx <- 1/delta * (incgam(alpha + delta, tau/beta)- incgam(alpha, tau/beta))

  #calculating c2
  e_1 <- log(beta) + (1/(incgam(alpha, tau/beta)))*(approx)

  return(e_1)
}

Q_function <- function(param, alpha.k, beta.k){
  alpha <- param[1]
  beta <- param[2]
  # sum((k-1)log(data)) + sum((k-1)*E_1)-
  #   sum(data/theta) + sum(E_2/theta) - m*k*log(theta)-
  #   m*log(k)
  tau <- 200
  
  q_eval <- -m*alpha*log(beta) - m*(lgamma(alpha))+ (alpha-1)*(sum(log(data)))- (1/beta)*sum(data)-
    W*(alpha-1)*C1(alpha.k, beta.k, tau) - (W/beta)*C2(alpha.k,beta.k,tau)
  
  return(q_eval)
  
}

init <- c(2,2)
params <- optim(par=init, alpha.k = init[1], beta.k = init[2], fn=Q_function,
                lower=c(0,0), 
                control=list(fnscale=-1))$par
params[1]*params[2]
mean(data)
