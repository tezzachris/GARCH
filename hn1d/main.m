
%Inputs needed:
%data: vector of log-returns
%rf: constant risk free rate

%Estimate parameters

start = []; %set empty if no idea about starting paramters
[a,ll, flag, fo] = fmincon_hn(data, rf, start);  

%Estimate st.errors

nw=0; %No newey west on scores
estimflag = 0; %not repeat estimation
[VCVrobust,A,~,scores,hess]=robustvcv('log_likelihood_hn',a,nw,data,estimflag,rf);
m = 1;
T = size(data,1);
VCV=hess^(-1)/(T-m);
se = sqrt(diag(VCV));
signif = ( abs( a ) ./ se ) >= tinv(0.95, T-length(a)) ;
