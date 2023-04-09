% backtest
% load VaRBacktestData
%     vbt = varbacktest(EquityIndex,...
%        [Normal95 Normal99 Historical95 Historical99 EWMA95 EWMA99],...
%        'PortfolioID','Equity',...
%        'VaRID',{'Normal95' 'Normal99' 'Historical95' 'Historical99' 'EWMA95' 'EWMA99'},...
%        'VaRLevel',[0.95 0.99 0.95 0.99 0.95 0.99]);



%1day-5day
%2004-2011
%from 1005 - 3005

%2012-2019
%from 3055 - 4555

%2019-end

rf = (0.13^2)/252; %risk free rate

primov= 1006; % 2004-01-06
ultimov= 4020; % 2015-12-24
subret = ret(primov:ultimov,1);

n = ultimov - primov + 1;
ninsample = ceil(0.5 * n);

j = 0;
nsim = 100000;
vars = zeros(n-ninsample,1); 
logliks = vars; fos = vars; exitflags=vars;
prevparam = [];
dayahead=2;
%h = waitbar(0,'Please wait...');

for i = ninsample : 1 : n-1
    j = j + 1 ; j
    %waitbar(i/(n-1),h)
    insample = subret(j:i,1);
    [param, loglik , exitflag, fo] = fminunc_chris2008(insample,rf,prevparam);
    if exitflag == 0 || fo > 0.1 %|| (param(7) + param(5)*param(6)^2)>=1 || param(2)*param(3)^2 + param(4) >=1
        prevparam = [];
        [param, loglik , exitflag, fo] = fminunc_chris2008(insample,rf,prevparam);
    else 
        prevparam = param;
    end
    %simuls = real( vecsimulate(2,param,rf,nsim) );
    simuls = vecsimulate_chris(dayahead,param,rf,nsim); %usavo 2 prima meglio 1?
    exitflags(j)=exitflag;
    logliks(j) = loglik;
    vars(j) = quantile(simuls(end,:), 0.05);
    fos(j) = fo;
    
end

%find(imag(simuls(2,:))~=0)
vbt=varbacktest(  subret(ninsample+1:n,1) , -vars,'VaRLevel',0.95);
kupiec=pof(vbt);


