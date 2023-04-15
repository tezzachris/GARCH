% backtest
% load VaRBacktestData
%     vbt = varbacktest(EquityIndex,...
%        [Normal95 Normal99 Historical95 Historical99 EWMA95 EWMA99],...
%        'PortfolioID','Equity',...
%        'VaRID',{'Normal95' 'Normal99' 'Historical95' 'Historical99' 'EWMA95' 'EWMA99'},...
%        'VaRLevel',[0.95 0.99 0.95 0.99 0.95 0.99]);



%1day-5day
%2004-2012
%1006 to 3020

%2012-2019
%5547 to 7308

%2019-end
%7308 to end


rf = (0.13^2)/252; %risk free rate

primov= 1006; % 2004-01-06
ultimov= 4020; % 2015-12-24

primov= 2006; % 2004-01-06
ultimov= 5020; % 2015-12-24
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
    [param, loglik , exitflag, fo] = fminunc_hn(insample,rf,prevparam);
    if exitflag == 0 || fo > 0.1
        prevparam = [];
        [param, loglik , exitflag, fo] = fminunc_hn(insample,rf,prevparam);
    else 
        prevparam = param;
    end
    simuls = vecsimulate_hn(dayahead,param,rf,nsim); %usavo 2 prima meglio 1?
    exitflags(j)=exitflag;
    logliks(j) = loglik;
    vars(j) = quantile(simuls(end,:), 0.05);
    fos(j) = fo;
    
end

vbt=varbacktest(  subret(ninsample+1:n,1) , -vars,'VaRLevel',0.95);
kupiec=pof(vbt); %proportion of failures
%runtests(vbt)

p=1-kupiec.VaRLevel; 
N=kupiec.Observations; %number of obs
x=kupiec.Failures; %failures (simulated VaR > or < empirical VaR)
%failures=sum( ret2(201:end,1) <  vars); deve essere uguale a x

% Confronta kupiec.LRatioPOF con seguente
LRpof = -2*log( ((1-p)^(N-x) * p^x) / ( (1-x/N)^(N-x) * (x/N)^x  ) );



