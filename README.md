## Affine GARCH Models (i.e. possess semi-analytical expression for characteristic function)

Each model has independent folder path that carries out the model: estimation, out-sample test, simulation, option pricing (via characteristic function. 
The models we study in the article [4] are:

1) Heston-Nandi in [1]
2) Engle-Lee in [2]
3) CJOW in [3]
4) F2 (with nested models)
Extra) CJOW Positive model

We briefly describe hereafter the content of each folder:

# Model Estimation 

For each model the estimation routine is composed of several functions, the steps in a nutshell are:

1) Creation of a grid of (simulated) starting values for each parameter
2) Among those in the grid, select the combination that maximises the log-likelihood
3) Use MATLAB optmizers i.e. fmincon and fminunc to find the MLEs

 Note: Ideally one should repeat the steps above multiple times to have some confidence in the results.

# Model Simulation 

# Option Pricing via Characteristic Function 

 ## References

 [1]
 [2]


