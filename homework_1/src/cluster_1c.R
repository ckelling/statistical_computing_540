###
### Problem 1
###

library(tictoc)
# Write computer code to simulate M = 100 draws from a multivariate normal
#     distribution with mean 0 (vector of 0s) and covariance matrix
#     as above for n = 1000.
#     Sigma_{ij} = exp(-|i-j|/phi),

#Creating 3 different algorithms


sim_alg1 <- function(n, mu, Sigma){ #cholesky factorization
  X <- matrix(rnorm(n*M), nrow = M, ncol=n)
  L <- chol(Sigma)
  mv_n <- matrix(mu, nrow = M, ncol =n) + X%*%L #even through mu is 0 in our case
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
### 1c) Now that you have successfully implemented this for n = 1000, repeat the exercise 
###      for n = 2,000; n = 3,000 and so on, going up to the largest value you are able 
###      to handle. Make plots of how the computational cost scales with increasing n.
###

n <- c(1:5*1000,11:14*500)
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
  save(sys_time2, file= "/storage/home/cek32/work/prob_1c_sys_clust.Rdata")
}




