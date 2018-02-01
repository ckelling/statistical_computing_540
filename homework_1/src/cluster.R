library(RSpectra)
setwd("/storage/home/cek32/work/")

m <- c(10000)
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
    d1_new <- lambda[3]-lambda[2]
    d2_new <- lambda[2]-lambda[1]
    d1[i] <- d1_new
    d2[i] <- d2_new
  }
  df <- cbind(df, d1,d2)
}

save(df, file = "/storage/home/cek32/work/df.Rdata")