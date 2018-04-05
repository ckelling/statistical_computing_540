prob_5 <- scan("http://personal.psu.edu/muh10/540/data/bulbsHW3.dat")
data <- as.vector(prob_5)

tau <- 200 #stopping time in number of days
m <- 300
obs <- length(data) #number of observed lifetimes
W <- m-obs #number of unobserved lifetimes

###
### Coding up Paper
###

C2 <- function(alpha, beta, tau){
  num <- beta*incgam(a = alpha+1,x= tau/beta)
  denom <- incgam(a= alpha,x= tau/beta)
  c_2 <- num/denom
  return(c_2)
}


C1 <- function(alpha, beta, tau){
  #approximating derivative
  delta <- 1e-7
  approx <- (1/delta) * (incgam(a = alpha + delta, x = tau/beta)- incgam(a = alpha,x = tau/beta))
  
  #calculating c2
  c_1 <- log(beta) + (1/(incgam(a = alpha, x = tau/beta)))*(approx)
  
  return(c_1)
}

# Testing finite difference method for this function
alpha <- 10
beta <- 15
tau <- 200

#numerical instability around 1e-10
#stick with delta around 1e-7
delta_vec <- seq(1e-10,1e-6, by=1e-10)
approx_vec <- NULL
for(i in 1:length(delta_vec)){
  approx_vec[i] <- (1/delta_vec[i])^2 * (incgam(a = alpha - delta_vec[i], x = tau/beta)- 2*incgam(a = alpha,x = tau/beta) + 
                                 incgam(a = alpha + delta_vec[i], x = tau/beta))
  #approx_vec[i] <- 1/delta_vec[i] * (incgam(a = alpha + delta_vec[i],x= tau/beta)- incgam(a =alpha,x= tau/beta))
}

i <- 1:1000
plot(delta_vec[-i], approx_vec[-i])
length(delta_vec)
delta_vec[500]



Q_function <- function(param, alpha.k, beta.k){
  alpha <- param[1]
  beta <- param[2]
  tau <- 200
  
  q_eval <- -m*alpha*log(beta) - m*(lgamma(alpha))+ (alpha-1)*(sum(log(data)))- (1/beta)*sum(data)+
    W*(alpha-1)*C1(alpha.k, beta.k, tau) - (W/beta)*C2(alpha.k,beta.k,tau)
  
  return(q_eval)
  
}


#test it out
init <- c(2,20)
params <- optim(par=init, alpha.k = init[1], beta.k = init[2], fn=Q_function,
                lower=c(1e-20,1e-20), 
                control=list(fnscale=-1))$par


f_em <- function(y, a.1, b.1, eps, maxit){
  #input:
  #    y: times
  #    eps: tolerance for the convergence of beta
  #    a: initial value
  #    b: initial value
  #    maxit: max number of iterations
  
  
  #just to calculate initial distance
  a <- Inf
  b <- Inf
  
  diff.a <- abs(a- a.1)
  diff.b <- abs(b- b.1)

  i <- 1
  EM.hist <- data.frame(i, diff.a, diff.b)
  param.hist <- matrix(c(a.1, b.1), nrow = 1)
  
  #   while we have not had convergence and are not at our maximum number of iterations....
  while((i <= maxit) & ((diff.a > eps) | (diff.b > eps))){
    i <- i+1
    # make them into old param
    a.2 <- a.1
    b.2 <- b.1
    
    #compute optimal
    init <- c(a.2,b.2)
    params <- optim(par=init, alpha.k = init[1], beta.k = init[2], fn=Q_function,
                    lower=c(1e-20,1e-20), 
                    control=list(fnscale=-1))$par
    ?optim
    
    #update param
    a.1 <- params[1]
    b.1 <- params[2]
    
    #now, calculate new distances (to assess for convergence)
    diff.a <- abs(a.2- a.1)
    diff.b <- abs(b.2- b.1)
    
    # iteration history
    EM.hist   <- rbind(EM.hist, c(i, diff.a, diff.b))
    param.hist <- rbind(param.hist, matrix(c(a.1, b.1), nrow = 1))
  }
  
  # prepare output
  #       beta.MLE: after convergence, what is beta_MLE
  #       iter: how long it took to converge
  #       NR.hist: history matrix as above, all information on convergence differences
  #       param.hist: history of betas through iterations
  output <- list()
  output$a.MLE  <- a.1
  output$b.MLE  <- b.1
  output$iter      <- i - 1
  output$EM.hist   <- EM.hist
  output$param.hist <- param.hist
  
  #calculate standard error
  
  return(output)
}

#initial values
a.init <- 10;    b.init <- 10
#a.init <- 100;  b.init <- 100
#a.init <- 160;  b.init <- 10
#a.init <- 50;   b.init <- 100

tic()
output_0 <- f_em(y, a.init, b.init, 1e-6, 700)
em_time <- toc()
#em_time <- em_time$toc-em_time$tic


output_0$iter
p_hist <- output_0$param.hist
em_hist <- output_0$EM.hist

param_df <- cbind(c(rep("alpha", nrow(output_0$param.hist)),rep("beta", nrow(output_0$param.hist))),
                 c(output_0$param.hist[,1],output_0$param.hist[,2] ))
colnames(param_df) <- c("parameter", "est")
param_df <- as.data.frame(param_df)
param_df$ind <- c(1:nrow(output_0$param.hist),1:nrow(output_0$param.hist))
param_df$est <- as.numeric(as.character(param_df$est))

ggplot(data=param_df, aes(x=ind,y=est, group=parameter, col = parameter)) +
  geom_line()+
  geom_point()+labs(title = "EM Estimation with Initial Values (2, 2)")



#MLE Values
output_0$a.MLE
output_0$b.MLE

#for calculating standard errors
alpha <- output_0$a.MLE
beta <- output_0$b.MLE

#Standard Error
I_T_b <- (m*alpha)/(beta^2)
I_T_a <- (m*trigamma(alpha))


exp_c <- (beta* incgam(a = alpha+1, x = tau/beta))/(incgam(a = alpha, x = tau/beta))

I_c_b <- (2*exp_c)/(beta^3) - (alpha/beta^2) + 
  (exp(-tau/beta)*tau^alpha)/(incgam(a=alpha, x =tau/beta)*beta^(alpha+2))*((tau/beta)-alpha-1)-
  ((exp(-tau/beta)*(tau^alpha)/(beta^(alpha+1)))/(incgam(a=alpha,x=tau/beta)))^2


delta <- 6e-7
first_approx <- (1/delta) * (incgam(a = alpha + delta, x = tau/beta)- incgam(a = alpha,x = tau/beta))
sec_approx <- (1/delta)^2 * (incgam(a = alpha - delta, x = tau/beta)- 2*incgam(a = alpha,x = tau/beta) + 
                               incgam(a = alpha + delta, x = tau/beta))


I_c_a <- (sec_approx)/(incgam(a=alpha, x =tau/beta)) - 
  ((first_approx)/(incgam(a=alpha, x=tau/beta)))^2

I_Y_a <- I_T_a - I_c_a
I_Y_b <- I_T_b - I_c_b

inv_I_Y_a <- (I_Y_a)^-1
inv_I_Y_b <- (I_Y_b)^-1

#Confidence Intervals
alpha + 1.96*inv_I_Y_a/sqrt(m)
alpha - 1.96*inv_I_Y_a/sqrt(m)

beta + 1.96*inv_I_Y_b/sqrt(m)
beta - 1.96*inv_I_Y_b/sqrt(m)
