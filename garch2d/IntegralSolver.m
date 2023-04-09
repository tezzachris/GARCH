% GetHesnanPrice.m
%function [cHN]=GetHesnanPrice(S,T,K,r,a,b,om,gamstar,ht1,resc);
% Parameter Values mostly from Heston and Nandi (2000) Table 1(a)
r=0; % risk free rate per day continuously compounded 
T=90; %days to maturity
a=1.32e-6; % GARCH MA parameter
b=0.589; % GARCH AR parameter
gam=421.39; % GARCH Asymmetry parameter
om=5.02e-6; % GARCH intercept
lam=.205; % risk premium
ht1=(0.15^2)/252; % realistic value for daily unconditional variance 
K=100;% strike price
S=100; % stock price
gamstar = gam + .5 + lam;

integr=integral(@(x) hnint(x,a,b,om,ht1,-.5,gamstar,T,r,K,S),0,+Inf,'ArrayValued',true); 
% Calculate Heston European call price
op_price=0.5*(S-K*exp(-r*T)) + K/pi*integr*exp(-r*T)


