

%needs function varbacktest
%three functions: fmincon, fminunc, simulate

function [vbt1,QL1,vbt5,QL5,logliks] = outsample_hn(primov, ultimov, ret, rf)

subret = ret(primov:ultimov,1);
n = ultimov - primov + 1;
ninsample = ceil(0.6 * n);
j = 0; %indicator for the loop
nsim = 100000; %as paper enzo/ballestra
VARs1 = zeros(n-ninsample,1);  %storing the value at risk
VARs2=VARs1; VARs3=VARs1; VARs4=VARs1; VARs5=VARs1; 
prevparam = [];
logliks = zeros(n-ninsample,1);

for i = ninsample : 1 : n-1
    j = j + 1 ; j
    insample = subret(j:i,1);
    [param, loglik , exitflag,fo] = fmincon_hn(insample,rf,prevparam);

    while exitflag <= 0 
        prevparam = [];
        [param, loglik , exitflag] = fmincon_hn(insample,rf,prevparam);
    end 
    logliks(j) = loglik;
    prevparam = param;
    simuls = vecsimulate_hn(6,param,rf,nsim); %some can be complex number as variance is not guarantee positive
    VARs1(j) = quantile( simuls(2,:), 0.05);
    VARs2(j) = quantile( simuls(3,:), 0.05);
    VARs3(j) = quantile( simuls(4,:), 0.05);
    VARs4(j) = quantile( simuls(5,:), 0.05);
    VARs5(j) = quantile( simuls(6,:), 0.05);

end
vbt1=varbacktest(  subret(ninsample+1:n,1) , -VARs1, 'VaRLevel',0.95);
vbt2=varbacktest(  subret(ninsample+1:n,1) , -VARs2, 'VaRLevel',0.95);
vbt3=varbacktest(  subret(ninsample+1:n,1) , -VARs3, 'VaRLevel',0.95);
vbt4=varbacktest(  subret(ninsample+1:n,1) , -VARs4, 'VaRLevel',0.95);
vbt5=varbacktest(  subret(ninsample+1:n,1) , -VARs5, 'VaRLevel',0.95);

yt = subret(ninsample+1:end);
QL1 =  ( (yt - VARs1) .* (0.05 - (yt < VARs1) ) );
QL5 =  ( (yt - VARs5) .* (0.05 - (yt < VARs5) ) );

end


