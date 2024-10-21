

%Log-returns
[a,ll] = fminunc_hn(ret, rf, []);

%Options
[a,ll] = fmincon_hn_options(ret, df, rf );


%Standarad errors
nw=0; %No newey west on scores
estimflag = 0;
[VCVrobust]=robustvcv('log_likelihood_hn',a,nw,ret,estimflag,rf);
se = sqrt(diag(VCVrobust)); 
sig = ( abs( a ) ./ se ) >= tinv(0.95, size(ret,1)-length(a)) ;
