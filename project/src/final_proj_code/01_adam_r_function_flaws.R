###
### Claire Kelling
### SODA 540
### Final Project
###
### Created 3/8/18 for final project code on Adam optimization
### Exploring the built-in R function that implements Adam optimization.
### Many flaws, listed at the bottom of this file.
### 

# Libraries
library(gradDescent)

# Data
## Learning and Build Model with ADAM, Example from R package (for reference)
## load R Package data
data(gradDescentRData)

# This is a dataset that collected experimentaly by Kennedy in 1954 to obtain the density value of
# CO2.  Parameter used in the experiment are temperature and pressure, which it can be a parameter
# to obtain compressibility factor value of CO2.

## get z-factor data
dataSet <- gradDescentRData$CompressilbilityFactor

## split dataset into training and testing data
split_dat <- splitData(dataSet)

## build model with ADAM
ADAMmodel <- ADAM(split_dat$dataTrain)

#show result, a vector matrix of theta (coefficient) for linear model
# (no other output other than the vector matrix of theta)
# (how to assess for convergence?)
# Also, this algorithm is not very stable in between training sets or in between runs.
print(ADAMmodel)

# Try fitting linear model to compare the results (pretty different than ADAM results)
# The coefficients do not appear to make a lot of sense in the linear model framework.
# This is because there is a pre-processing to the data, where there is normalization of the data.
tr_dat <- split_dat$dataTrain
colnames(tr_dat) <- c("temp", "pressure", "comp_fac")
ln_mod <- lm(comp_fac ~ temp + pressure, data = tr_dat)
summary(ln_mod)$coefficients


##
## Another example of model fit
##
dataSet <- gradDescentRData$CompressilbilityFactor
## train dataset
modelObject <- gradDescentR.learn(dataSet)

temp <- c(273.1, 353.1, 363.1)
pres <- c(24.675, 24.675, 24.675)
conf <- c(0.8066773, 0.9235751, 0.9325948)
zfac <- data.frame(temp, pres, conf)
## predict
prediction_data <- predict(modelObject, zfac)

##
## Example of pre-processing the data
##
square_feet <- c(1400,1600,1700,1875,1100,1550,2350,2450,1425,1700)
price <- c(245,312,279,308,199,219,405,324,319,255)

house_price <- data.frame(square_feet, price)

## Preprocessing stage
# This function is degraded- no longer workds
#house_price.data <- gradDescent::gradDescent.preprocessing(
#  house_price,
#  normalizationMethod="variance"
#)

house_price.data2 <- normalize(house_price)


## Model building stage
GD.model <- gradDescent.learn(
  house_price.data,
  methodType="GD" 
)

#
# Test pre-processing on the other dataset
#
# This function is degraded- no longer works
#split_norm <- gradDescent.preprocessing(
#  split_dat$dataTrain,
#  normalizationMethod="variance"
#)

split_norm <- normalize(split_dat$dataTrain)

ADAMmodel2 <- ADAM(split_norm)

#compare output of ADAM with and without pre-processing (they are similar, though unstable)
ADAMmodel
ADAMmodel2

#Now, we will compare to the coefficients returned by the linear model when the 
# data is normalized. We find they are somewhat similar.
colnames(split_norm) <- c("temp", "pressure", "comp_fac")
ln_mod <- lm(comp_fac ~ temp + pressure, data = split_norm)
summary(ln_mod)$coefficients


###
### Takeaways
###

# 1.) The R function ADAM returns "parameters for a linear model." However, these aren't close to
# the parameters returned by fitting a linear model. ADAM rescales the data, and returns the coefficients
# for the re-scaled data. 

# 2.) This function, for this data, does not return stable values. if we rerun the algorithm many 
# times, we get many different estimates for the coefficients. 

# 3.) The function also does not return anything that would allow me to assess for convergence, 
# such as iteration history.

# 4.) There are no inputs for this funtion for beta_1, beta_2 or epsilon. There are no listed
# default values in the documentation of the function, although they suggest a couple in the paper
# (but only for machine learning problems).