###
### Claire Kelling
### STAT 540
### Homework #3
###
### Created 3/27/2018 for the third assigment, due 4/5/2018
### The exercises mainly focus on ??
### 

###
### Problem 1
###

# Find the MLE for a logistic regression model.
# Write a Newton-Raphon algorithm to find the MLE (\hat{\beta_1}; \hat{\beta_2}). 
# The data for this problem are available here: http://personal.psu.edu/muh10/540/data/logReg.dat

#clear workspace
rm(list=ls())

#load data
y <- fread("http://personal.psu.edu/muh10/540/data/logReg.dat")

# (a) Provide details about how you set up the algorithm, including any
#  analytical quantities you derived. You should also include the stop-
#  ping criteria, how you obtained initial values, how you tuned sn (step
#  size/learning rate) and other relevant information, e.g. did you run
#  the algorithm using multiple initial values?


# (b) Provide your pseudocode for the main algorithm. You can use no-
#  tation defined in class (and specified above) in order to make the
#  pseudocode succinct and easy to read

# (c) Provide the MLE as well as estimate of the standard error, along
#  with asymptotic 95% confidence intervals.


# (d) Report the total computing time taken by your algorithm (e.g. you
#  could use system.time). If you have 2-3 different versions of the
#  algorithm, report timings for each version.


# (e) Provide the computational cost of the algorithm (flops) per iteration
#  as a function of n. Also provide the number of iterations it took before
#  the stopping criteria was met. Of course, the number of iterations will
#  vary by initial conditions so you should provide at least 2-3 different
#  counts depending on where you started.

# (f) Summarize in 1-2 sentences the sensitivity of the algorithm to the
#  initial value, that is, does the algorithm work well regardless of where
#  you start or does it need to be in a neighborhood of a certain value?
#    (This is not required, but because this is a simple 2-D optimization
#     problem, you could also provide the 2-D log-likelihood surface to
#     obtain insights.)

###
### Problem 2
###


# Implement a gradient descent algorithm (work with the negative log like-
# lihood). Answer (a)-(f) for this algorithm.

###
### Problem 3
###

# Now implement a stochastic gradient descent algorithm and repeat (a)-(f).

###
### Problem 4
###

#  Compare the three algorithms above: Newton-Raphson, gradient descent,
# and stochastic gradient descent. Provide a 1-2 sentence summary of which
# algorithm you would recommend and why. Then provide a more detailed
# comparison, for example how stable the algorithms are with respect to
# initial values, how sensitive they are to choice of sn, a comparison of the
# total computational cost, and the computational cost per iteration.


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

times <- fread("http://personal.psu.edu/muh10/540/data/bulbsHW3.dat")