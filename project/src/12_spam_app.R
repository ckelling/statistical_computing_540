###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/15/18 for final project use of spam dataset
### 

#https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.DOCUMENTATION
#https://archive.ics.uci.edu/ml/datasets/spambase
#https://www.kaggle.com/c/cs-189-logistic-regression2/data
spam_dat <- fread("https://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data")

source(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/src/00_opt_functions.R")

data <- spam_dat[,-(5:57)]
#data <- spam_dat

#try testing glm
mod <- glm(V58 ~.,family=binomial(link='logit'),data=spam_dat)
summary(mod)

beta1 <- 0.9
beta2 <- 0.999
eps1 <- 0.0000001 #for the algorithm
eps2 <- 1e-6 #convergence criteria
maxit <- 1000
step <- 1#step size of 1 helps speed

data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]

init <- rep(1, ncol(X))

sgd_spam <- sgd_opt(data, step, eps2, init, maxit)
sgd_spam <- nadam_bc(data, step, beta1, beta2, eps1, eps2, init, maxit)

sgd_spam$iter
sgd_spam$backtrack.iter

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

# [1,]  0.2277371
# [2,] -0.1290784
# [3,]  0.2052651
# [4,]  3.4976576