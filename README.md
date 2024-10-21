## Affine GARCH Models (i.e. possess semi-analytical expression for characteristic function)

Each model has independent folder path that carries out the model: estimation, out-sample test, simulation, option pricing (via characteristic function. 
The models we analyse in [4] are: 

1) Heston-Nandi in [1] 
2) Engle-Lee in [2] 
3) CJOW in [3]
4) GARCH-CPC in [4]
5) GARCH2F in [5] (with nested models) 

We briefly describe hereafter the content of each folder: 

# Model Estimation  

For each model the estimation routine is composed of several functions, the steps in a nutshell are: 

1) Gridsearch for optimal starting param values
2) MATLAB optmizers i.e. fmincon and fminunc to find the MLEs.

 Note: Ideally one should repeat the steps above multiple times to have some confidence in the results. 
 Note2: In my experience, fminunc tends to be more accurate but requires additional care for enforcing parameter constraints 

# Model Simulation  

The model simulation requires an initial value for the volatility and for the return process (optional).

# Option Pricing via Characteristic function

The option price of a European call/put on can be computed using the formula in [1] at Equation (11), that leverages the inversion formula of Gil-Pelaez.

# Data

The risk free rate is freely available for download at: https://fred.stlouisfed.org/series/TB3MS \
The S&P500 data (returns and option data) was downloaded via Refinitiv Eikon Datastream using a University Licence

 ## References

  [1]  S.L. Heston and S. Nandi. A Closed-Form GARCH Option Valuation Model. 2000 \
  [2]  R.F. Engle and G. Lee. A Long-Run and Short-Run Component Model of Stock Return Volatility. 1999 \
  [3]  P. Christoffersen, K. Jacobs, C. Ornthanalai, and Y. Wang. Option Valuation with Long-run and
Short-run Volatility Components. 2008  \
  [4] GARCH-CPC: https://arxiv.org/html/2410.14513v1
  [5] GARCH2F: https://arxiv.org/abs/2410.14585


