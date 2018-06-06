This is the source code file for my final project for STAT 540: Statistical Computing.

Each file represents a step in the analysis. I believe that it is important to code this way 
in order to increase reproducibility at each step of the code. Also, it is easier to document
each individual step of the analysis this way.

Please find brief descriptions of each file below. 

00_opt_functions.R: Functions that will later be used to analyze the spambase dataset. This
	includes all of the optimization algorithms as well as the objective and score functions.

01_adam_r_functions_flaws.R: Exploring the built-in R function that implements Adam optimization.
	I have discovered many flaws of the existing implementation of Adam, documented at the 
	end of the code file.

02_adam_opt_lin_reg.R: Coding the adam algorithm for linear regression.

03_adam_opt_logreg.R: Coding the adam algorithm for logistic regression, as an exploration. This
	does not include detailed storage mechanisms for later analysis.

04_sgd.R: A detailed implementation of SGD, which then stores the information for later analysis.

05_sgdm.R: A detailed implementation of SGD with momentum, which then stores the information 
	for later analysis.

06_adam.R: A detailed implementation of Adam, which then stores the information for later analysis.

07_nadam.R: A detailed implementation of Nadam, which then stores the information for later analysis.

08_nobc_adam.R: A detailed implementation of Adam without the bias correction steps, which then 
	stores the information for later analysis.

09_nobc_nadam.R: A detailed implementation of Nadam without the bias correction steps, which then 
	stores the information for later analysis.

10_nag.R: A detailed implementation of NAG, which then stores the information for later analysis.

11_combined_analysis.R: Combines the results of the preceding 7 algorithms and plots the convergence
	of both the parameters and the results.

12_spam_app.R: Exploration of the spam dataset.

13_full_spam_analysis.R: Completing the optimization steps for the spam dataset, using the 00_opt_functions.R
	file for the optimization functions.

14_spam_plot.R: Plotting the results from the previous file, to see convergence of both the objective
	function as well as the parameters.