# Estimate some GARCH models

For each model the estimation routine is composed of several functions, the steps in a nutshell are:

1) Creation of a grid of (simulated) starting values for each parameter
2) Among those in the grid, select the combination that maximises the log-likelihood
3) Use MATLAB optmizers i.e. fmincon and fminunc to find the MLEs

 Note: Ideally one should repeat the steps above multiple times to have some confidence in the results.


