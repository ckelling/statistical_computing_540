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
             mean(nag_store$df$time),mean(adam_store$df$iter),
             mean(nadam_store$df$iter),mean(adam_store_nobc$df$iter),
             mean(nadam_store_nobc$df$iter))

citer_df <- c(mean(sgd_store$df$c_iter),mean(sgdm_store$df$c_iter),
             mean(nag_store$df$time),mean(adam_store$df$c_iter),
             mean(nadam_store$df$c_iter),mean(adam_store_nobc$df$c_iter),
             mean(nadam_store_nobc$df$c_iter))

titer_df <- c(mean(sgd_store$df$tot_iter),mean(sgdm_store$df$tot_iter),
              mean(nag_store$df$time),mean(adam_store$df$tot_iter),
              mean(nadam_store$df$tot_iter),mean(adam_store_nobc$df$tot_iter),
              mean(nadam_store_nobc$df$tot_iter))

full_comp <- rbind(time_df, iter_df, citer_df, titer_df)
colnames(full_comp) <- c("sgd","sgd_m","nag","adam", "nadam", "adam_nobc", "nadam_nobc")
rownames(full_comp) <- c("time", "iter", "citer", "titer")

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

plot_dat$coeff2 <- factor(plot_dat$coeff, labels = c("theta[1]", "theta[2]"))

ggplot(data=plot_dat, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)

#Earlier Cutoff for graph:
half <- nrow(sgd_df)/2+1
n <- 350
plot_dat2 <- rbind(sgd_df[c(1:n, half:(half+n-1)),],sgdm_df[c(1:n, half:(half+n-1)),],
                   nag_df[c(1:n, half:(half+n-1)),],adam_df[c(1:n, half:(half+n-1)),], 
                   nadam_df[c(1:n, half:(half+n-1)),], adam_df_nobc[c(1:n, half:(half+n-1)),],
                   nadam_df_nobc[c(1:n, half:(half+n-1)),])
plot_dat2$coeff2 <- factor(plot_dat2$coeff, labels = c("theta[1]", "theta[2]"))
colnames(plot_dat2)[1] <- "Algorithm"

ggplot(data=plot_dat2, aes(x=ind,y=est, group=Algorithm, col = Algorithm)) +
  geom_line(size = 1)+
  #geom_point()+
  labs(title = "Convergence Comparison")+
  #facet_wrap(~ coeff2, ncol = 3, labeller = label_parsed)
  facet_grid(.~coeff2, labeller = label_parsed)
