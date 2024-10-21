
% Code for estimating the GARCH2F model.

% Log-returns

% Requires: vector of log-returns(ret) and risk-free rate (rf)

initial = []
estim_flag = 0 %flag = 0 fmincon, = 1 fminunc 
[a, ll] = fmincon_full(ret, rf, initial); [val, v, h] = loglik_full(a,ret,estim_flag,rf);
[a, ll] = fminunc_full(ret, rf, initial);

% Options

%Requires: table of option prices, strikes, days to maturity, flag call/put, underlying

[a,ll] = fmincon_hn_options(ret, df, rf );


% Standard errors

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Date: 9/1/2005

nw=0; %No newey west on scores 
estimflag = 0; 
[VCVrobust]=robustvcv('loglik_full',a,nw,ret,estimflag,rf); 
se = sqrt(diag(VCVrobust)); %standard errors 
signif = 0.95 
sig = ( abs( a ) ./ se ) >= tinv(signif, length(ret)-length(a)) ;

