###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/10/18 for final project comparison of all algorithms considered.
### 

#sgd
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_sgd_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/sgd_df.Rdata")
#sgd with momentum
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_sgdm_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/sgdm_df.Rdata")
#nag
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nag_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nag_df.Rdata")
#adam
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nadam_sim2.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/adam_df.Rdata")
#nadam
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nadam_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_df.Rdata")
#adam, no bias correction
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_adam_nobc_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/adam_df_nobc.Rdata")
#nadam, no bias correction
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/full_nadam_nobc_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nadam_df_nobc.Rdata")


## Comparison in terms of Computation Time:
# Have to convert to dataframe first
sgd_store$df <- as.data.frame(sgd_store$df)
sgdm_store$df <- as.data.frame(sgdm_store$df)
nag_store$df <- as.data.frame(nag_store$df)
adam_store$df <- as.data.frame(adam_store$df)
nadam_store$df <- as.data.frame(nadam_store$df)
adam_store_nobc$df <- as.data.frame(adam_store_nobc$df)
nadam_store_nobc$df <- as.data.frame(nadam_store_nobc$df)

time_df <- c(mean(sgd_store$df$time),mean(sgdm_store$df$time),
             mean(nag_store$df$time),mean(adam_store$df$time),
             mean(nadam_store$df$time),mean(adam_store_nobc$df$time),
             mean(nadam_store_nobc$df$time))

## Comparison in terms of Number of iterations:
iter_df <- c(mean(sgd_store$df$iter),mean(sgdm_store$df$iter),
             mean(nag_store$df$iter),mean(adam_store$df$iter),
             mean(nadam_store$df$iter),mean(adam_store_nobc$df$iter),
             mean(nadam_store_nobc$df$iter))

citer_df <- c(mean(sgd_store$df$c_iter),mean(sgdm_store$df$c_iter),
             mean(nag_store$df$c_iter),mean(adam_store$df$c_iter),
             mean(nadam_store$df$c_iter),mean(adam_store_nobc$df$c_iter),
             mean(nadam_store_nobc$df$c_iter))

titer_df <- c(mean(sgd_store$df$tot_iter),mean(sgdm_store$df$tot_iter),
              mean(nag_store$df$c_iter),mean(adam_store$df$tot_iter),
              mean(nadam_store$df$tot_iter),mean(adam_store_nobc$df$tot_iter),
              mean(nadam_store_nobc$df$tot_iter))

full_comp <- rbind(time_df, iter_df, citer_df, titer_df)
colnames(full_comp) <- c("sgd","sgd_m","nag","adam", "nadam", "adam_nobc", "nadam_nobc")
rownames(full_comp) <- c("time", "iter", "citer", "titer")

full_comp_report <- full_comp[c(1,2),]
xtable(full_comp_report)

#Quick fix on the name of the nadam_nobc algorithm
nadam_df_nobc$algo <- rep("nadam_nobc", nrow(nadam_df_nobc))

## Comparison in terms of Computation Time per iteration:
# iter_t_df <- c(mean(sgd_store$df$tot_iter/), mean(adam_store$df$tot_iter),
#               mean(nadam_store$df$tot_iter),mean(adam_store_nobc$df$tot_iter),
#               mean(nadam_store_nobc$df$tot_iter))

## Comparison in terms of Convergence:
# Need to combine all of the dataframes for plotting
plot_dat <- rbind(sgd_df, sgdm_df,nag_df, adam_df, nadam_df, adam_df_nobc, nadam_df_nobc)
colnames(plot_dat)[1] <- "Algorithm"

plot_dat$coeff2 <- factor(plot_dat$coeff, labels = c("beta[1]", "beta[2]"))

ggplot(data=plot_dat, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)

#Earlier Cutoff for graph:
half <- nrow(sgd_df)/2+1
n <- 250
plot_dat2 <- rbind(sgd_df[c(1:n, half:(half+n-1)),],sgdm_df[c(1:n, half:(half+n-1)),],
                   nag_df[c(1:n, half:(half+n-1)),],adam_df[c(1:n, half:(half+n-1)),], 
                   nadam_df[c(1:n, half:(half+n-1)),], adam_df_nobc[c(1:n, half:(half+n-1)),],
                   nadam_df_nobc[c(1:n, half:(half+n-1)),])
plot_dat2$coeff2 <- factor(plot_dat2$coeff, labels = c("beta[1]", "beta[2]"))
colnames(plot_dat2)[1] <- "Algorithm"

ggplot(data=plot_dat2, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison", x = "iterations")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)


#also plotting objective function
##Plotting the objective function over the averaged iteration
## we need to create this structure for all iterations of the algorithm
## we can just load the data below (instead of running below)
obj_fun = function(theta){
  x.b = as.vector(X%*%theta)
  return( -sum(y*x.b - log( 1 + exp(x.b))))  
}
prob1 <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")
data <- prob1
init <- c(10,10)
data<- as.data.frame(data)
X <- data[,1:(ncol(data)-1)]
X <- as.matrix(X)
y <- data[,ncol(data)]

## Start here#######################################
store <- nadam_store_nobc
obj_fun_st <- matrix(NA,nrow= nrow(store$theta.1.hist), ncol=ncol(store$theta.1.hist))
#information on convergence of objective function
for(j in 1:ncol(store$theta.1.hist)){
  for(i in 1:nrow(store$theta.1.hist)){
    #i <- 1
    #j<-1
    obj_fun_st[i,j] <- obj_fun(c(store$theta.1.hist[i,j],store$theta.2.hist[i,j]))
  }
}

obj_df <- NULL
obj_df$obj <- rowMeans(obj_fun_st)
obj_df <- as.data.frame(obj_df)
obj_df <- cbind(c(rep("nadam_nobc", nrow(obj_df))),
                c(obj_df$obj))
obj_df <- as.data.frame(obj_df)
colnames(obj_df) <- c("algo", "est")
obj_df$ind <- c(1:nrow(obj_df))
obj_df$est <- -as.numeric(as.character(obj_df$est))

nbc_obj2 <- obj_df


#full combination
full_obj2 <- rbind(s_obj2, sm_obj2, nag_obj2,
                  a_obj2, n_obj2,
                  abc_obj2, nbc_obj2)

#full_obj <- rbind(s_obj, a_obj, n_obj, abc_obj, nbc_obj)

#shorter
n <- 200
full_obj2 <- rbind(s_obj2[1:n,], sm_obj2[1:n,], nag_obj2[1:n,],
                   a_obj2[1:n,], n_obj2[1:n,],
                   abc_obj2[1:n,], nbc_obj2[1:n,])

#plotting the objective function convergence
ggplot(data=full_obj2, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence of Objective Function") #save 700 x 400
