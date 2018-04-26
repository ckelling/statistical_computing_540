###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/15/18 for final project use of spam dataset
### In this file, we just document the dataset, and do some initial exploration
### 

library(data.table)
library(ggplot2)

#https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.DOCUMENTATION
#https://archive.ics.uci.edu/ml/datasets/spambase
#https://www.kaggle.com/c/cs-189-logistic-regression2/data
spam_dat <- fread("https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data")

source(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/src/00_opt_functions.R")


#We explore different subsets of the data to use for modeling
#data <- spam_dat[,-c(4, 8, 9, 10, 11, 13:57)] #[,-c(4, 7:57)]
#data <- spam_dat[,c(3, 5, 16, 19, 52, 50, 55, 56, 57, 58)]
#data <- spam_dat[,-c(4, 7:57)]
#data <- spam_dat[,c(1,2,3,5,6,58)] #also 8, 9,10,11,12,13,14
data <- spam_dat[,c(1,2,3,5,6,8,9,10,11,12,13,14,15,58)] #now we have 13 covariates
mod <- glm(V58 ~.,family=binomial(link='logit'),data=data)
summary(mod)

#assess for multicolinearity
#pairs(data)
#data <- spam_dat

#try testing glm
#need to also include an intercept in our full analysis
mod <- glm(V58 ~.,family=binomial(link='logit'),data=data)
summary(mod)

#also testing glm for the other problem, with only two covariates
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
mod2 <- glm(y ~.,family=binomial(link='logit'),data=prob1)
summary(mod2)

#values for the algorithm
beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 10000
step <- 1#step size of 1 helps speed

#loading the data setup for the algorithm
data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- cbind(rep(1, nrow(data)), X) #also include an intercept
X <- as.matrix(X)
y <- data[,ncol(data)]

#initial values for the paramters
init <- rep(1, ncol(X))

#testing out convergence for the analysis
sgd_spam <- sgd_opt(data, step, eps2, init, maxit)
#sgd_spam <- nadam_bc(data, step, beta1, beta2, eps1, eps2, init, maxit)

#testing for convergence criteria
sgd_spam$iter
sgd_spam$backtrack.iter

#plotting the converged values
theta_hist <- sgd_spam$theta.hist

theta_df <- cbind(c(rep("theta_1", nrow(sgd_spam$theta.hist)),rep("theta_2", nrow(sgd_spam$theta.hist))),
                  c(sgd_spam$theta.hist[,1],sgd_spam$theta.hist[,2] ))
colnames(theta_df) <- c("coeff", "est")
theta_df <- as.data.frame(theta_df)
theta_df$ind <- c(1:nrow(sgd_spam$theta.hist),1:nrow(sgd_spam$theta.hist))
theta_df$est <- as.numeric(as.character(theta_df$est))

ggplot(data=theta_df, aes(x=ind,y=est, group=coeff, col = coeff)) +
  geom_line()+
  geom_point()+labs(title = "Stochastic Gradient Descent")
sgd_spam$theta.final

summary(mod)

#creating a storage mechanism for the objective function
obj_hist <- rep(NA, nrow(sgd_spam$theta.hist))

for(i in 1:nrow(sgd_spam$theta.hist)){
  #this plots the likelihood of the various theta iterations
  obj_hist[i] <- -obj_fun(sgd_spam$theta.hist[i,])
}

#plot the convergence of the objective function
plot(obj_hist)
