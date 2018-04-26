###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/20/18 for final project comparison of all algorithms considered.
### 

##Plotting the objective function over the averaged iteration
obj_df <- NULL
obj_df$obj <- rowMeans(store$obj_st)
obj_df <- as.data.frame(obj_df)
obj_df <- cbind(c(rep("nadam_nobc", nrow(obj_df))),
                     c(obj_df$obj))
obj_df <- as.data.frame(obj_df)
colnames(obj_df) <- c("algo", "est")
obj_df$ind <- c(1:nrow(obj_df))
obj_df$est <- -as.numeric(as.character(obj_df$est))

#s_obj <- obj_df
nbc_obj <- obj_df
#save(s_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/s_obj.Rdata")
#save(sm_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/sm_obj.Rdata")
#save(nag_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nag_obj.Rdata")
#save(a_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/a_obj.Rdata")
#save(n_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/n_obj.Rdata")
#save(abc_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/abc_obj.Rdata")
#save(nbc_obj,file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/nbc_obj.Rdata")

#quick correction
s_obj$algo <- rep("sgd", nrow(s_obj))

#full combination
full_obj <- rbind(s_obj, sm_obj, nag_obj,
                  a_obj, n_obj,
                  abc_obj, nbc_obj)

ggplot(data=full_obj, aes(x=ind,y=est, group=algo, col = algo)) +
  geom_line()+
  #geom_point()+
  labs(title = "Convergence of Objective Function")



###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 4/10/18 for final project comparison of all algorithms considered.
### 

#sgd
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_sgd_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_sgd_df.Rdata")
#sgd with momentum
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_sgdm_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_sgdm_df.Rdata")
#nag
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nag_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nag_df.Rdata")
#adam
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_sim2.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_adam_df.Rdata")
#nadam
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nadam_df.Rdata")
#adam, no bias correction
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_adam_nobc_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_adam_df_nobc.Rdata")
#nadam, no bias correction
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_full_nadam_nobc_sim.Rdata")
load(file = "C:/Users/ckell/OneDrive/Penn State/2017-2018/01_Spring/540/statistical_computing_540/project/data/spam_nadam_df_nobc.Rdata")


## Comparison in terms of Computation Time:
# Have to convert to dataframe first
s_store$df <- as.data.frame(s_store$df)
sm_store$df <- as.data.frame(sm_store$df)
nag_store$df <- as.data.frame(nag_store$df)
a_store$df <- as.data.frame(a_store$df)
n_store$df <- as.data.frame(n_store$df)
abc_store$df <- as.data.frame(abc_store$df)
nbc_store$df <- as.data.frame(nbc_store$df)

time_df <- c(mean(s_store$df$time),mean(sm_store$df$time),
             mean(nag_store$df$time),mean(a_store$df$time),
             mean(n_store$df$time),mean(abc_store$df$time),
             mean(nbc_store$df$time))

## Comparison in terms of Number of iterations:
iter_df <- c(mean(s_store$df$iter),mean(sm_store$df$iter),
             mean(nag_store$df$iter),mean(a_store$df$iter),
             mean(n_store$df$iter),mean(abc_store$df$iter),
             mean(nbc_store$df$iter))

citer_df <- c(mean(s_store$df$c_iter),mean(sm_store$df$c_iter),
              mean(nag_store$df$c_iter),mean(a_store$df$c_iter),
              mean(n_store$df$c_iter),mean(abc_store$df$c_iter),
              mean(nbc_store$df$c_iter))

titer_df <- c(mean(s_store$df$tot_iter),mean(sm_store$df$tot_iter),
              mean(nag_store$df$c_iter),mean(a_store$df$tot_iter),
              mean(n_store$df$tot_iter),mean(abc_store$df$tot_iter),
              mean(nbc_store$df$tot_iter))

full_comp <- rbind(time_df, iter_df, citer_df, titer_df)
colnames(full_comp) <- c("sgd","sgd_m","nag","adam", "nadam", "adam_nobc", "nadam_nobc")
rownames(full_comp) <- c("time", "iter", "citer", "titer")

full_comp_report <- full_comp[c(1,2),]
xtable(full_comp_report)

#Quick fix on the name of the nadam_nobc algorithm
#nadam_df_nobc$algo <- rep("nadam_nobc", nrow(nadam_df_nobc))

## Comparison in terms of Computation Time per iteration:
# iter_t_df <- c(mean(s_store$df$tot_iter/), mean(a_store$df$tot_iter),
#               mean(n_store$df$tot_iter),mean(abc_store$df$tot_iter),
#               mean(nbc_store$df$tot_iter))

## Comparison in terms of Convergence:
# Need to combine all of the dataframes for plotting
plot_dat <- rbind(s_storage_df, sm_storage_df,nag_storage_df, a_storage_df, n_storage_df, abc_storage_df, nbc_storage_df)
colnames(plot_dat)[1] <- "Algorithm"

plot_dat$coeff2 <- factor(plot_dat$coeff, labels = c("theta[1]", "theta[2]"))

ggplot(data=plot_dat, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)

#Earlier Cutoff for graph:
half <- nrow(nbc_storage_df)/2+1
n <- 600
m <- 600
plot_dat2 <- rbind(s_storage_df[c(1:m, half:(half+n-1)),],sm_storage_df[c(1:m, half:(half+n-1)),],
                   nag_storage_df[c(1:m, half:(half+n-1)),],a_storage_df[c(1:m, half:(half+n-1)),], 
                   n_storage_df[c(1:m, half:(half+n-1)),], abc_storage_df[c(1:m, half:(half+n-1)),],
                   nbc_storage_df[c(1:m, half:(half+n-1)),])
plot_dat2$coeff2 <- factor(plot_dat2$coeff, labels = c("theta[1]", "theta[2]"))
colnames(plot_dat2)[1] <- "Algorithm"

ggplot(data=plot_dat2, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison", x = "iterations")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)
