## GARCH2F Model 

# Model Estimation  

The estimation is composed of several functions: 

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
 GARCH2F: https://arxiv.org/abs/2410.14585


