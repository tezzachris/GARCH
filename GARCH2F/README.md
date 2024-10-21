
Code for estimating the GARCH2F model.

Example of use: 

%Requires: vector of log-returns(ret) and risk-free rate (rf)  

initial = []  
estim_flag = 0 %flag = 0 fmincon, = 1 fminunc - for parameter conversion  
[a, ll] = fmincon_full(ret, rf, initial); 
[val, v, h] = loglik_full(a,ret,estim_flag,rf); 

nw=0; %No newey west on scores 
estimflag = 0; 
[VCVrobust]=robustvcv('loglik_full',a,nw,ret,estimflag,rf); 
se = sqrt(diag(VCVrobust));  %standard errors 
signif = 0.95 
sig = ( abs( a ) ./ se ) >= tinv(signif, length(ret)-length(a)) ; 

