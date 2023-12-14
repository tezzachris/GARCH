# How to estimate *efficiently* GARCH models

# Aim: Find MLEs by the following numerical routine 

For each model the routine is composed of several functions, the steps in a nutshell are:
1) Creation of a grid of (simulated) starting values for each model parameter 
2) Using the "best" initial parameters, rely MATLAB function optmizers i.e. fmincon and fminunc to find MLEs

 Note: the MATLAB optimizers are "black-box" machinery so ideally one should repeat the steps above multiple times to have some confidence in the results.


