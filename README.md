## GARCH2F Model 

Garch model which includes two innovation factors and two volatility components

# Model Estimation  

The estimation of parameters can be tricky given the model complexity. Currently it is implemented in MATLAB and thus it relies its (non-linear) optimization functions, fmincon and fminunc. To facilitate parameter search I have created a starting parameter function to perform a standard gridsearch for optimal initial parameter values.

Note: Ideally one should repeat the steps above multiple times to have some confidence in the results. 
In my experience, fminunc tends to be more accurate but requires additional care for enforcing parameter constraints. 

# Model Simulation  

The model simulation requires an initial value for the volatility and for the return process (optional).

# Option Pricing via Characteristic function

The option price of a European call/put can be computed using the formula in [1] at Equation (11), that leverages the inversion formula of Gil-Pelaez.

# Data

The risk free rate is freely available for download at: https://fred.stlouisfed.org/series/TB3MS \
The S&P500 data (returns and option data) was downloaded via Refinitiv Eikon Datastream using a University Licence

 ## References
 GARCH2F: https://arxiv.org/abs/2410.14585


