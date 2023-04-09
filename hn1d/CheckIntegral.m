% GetHesnanPrice.m
% function [cHN]=GetHesnanPrice(S,T,K,r,a,b,om,gamstar,ht1,resc);
% Parameter Values mostly from Heston and Nandi (2000) Table 1(a)
r=0; % risk free rate per day continuously compounded T=90; % days (periods) to maturity
T=90;
a=1.32e-6; % GARCH MA parameter
b=0.589; % GARCH AR parameter
gam=421.39; % GARCH Asymmetry parameter
om=5.02e-6; % GARCH intercept
lam=.205; % risk premium
ht1=(0.15^2)/252; % realistic value for daily unconditional variance 
K=100;% strike price
S=100; % stock price
gamstar = gam - .5;




L=2e6;
logst=zeros(L,1);

for j = 1:L
    prova=simul_hn_1d(T);
    logst(j)=prova(end);
end

mean( max( exp(logst) - K ,0) )


