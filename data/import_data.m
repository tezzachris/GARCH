

%SPX500 Daily levels Adjuste Close Price
%2000-01-03
%2022-11-08
rf = (0.13^2)/252; %risk free rate
rawdata=importdata('sp500.csv');
lvl=rawdata.data;
ret = log( lvl(2:end) ./ lvl(1:end-1) )  ;
date=rawdata.textdata(2:end,2); %escludi nome colonna
%plot(datetime(date,ret))



djilvl=importdata('dji.csv');
willvl=importdata('wilsh.csv');

splvl=splvl.data; djilvl = djilvl.data; willvl=willvl.data; 

data = log( splvl(2:end) ./ splvl(1:end-1) )  ;
data1= log( djilvl(2:end) ./ djilvl(1:end-1) )  ;
data2= log( willvl(2:end) ./ willvl(1:end-1) ) ;
data2 = data2(2500:end);





%simulate a path with ht conditional variance
%simul_data = tarch_simulate(T,initial,1,0,1,[],2);

%smart function to detect optimal starting values

startingvals=[];
%data is zero mean residuals vector
fdata = data.^2;
fIdata = fdata.*(data<0);
p=1; o = 0 ; q=1;
T=length(data);
error_type='STUDENTST';
tarch_type= 2;

%params=tarch_starting_values(startingvals,data,fdata,fIdata,p,o,q,T,error_type,tarch_type);

%tparams=tarch_transform(params, p , o , q, 1);
%iparams=tarch_itransform(tparams, p , o , q, 1);



[parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = tarch(data, 1, o, q, error_type, tarch_type, startingvals, []);

back_cast= (theta(1))/(1-theta(2)-theta(3));
estim_flag = 1;

tarch_likelihood(theta, data, fdata, fIdata, p , o ,q, error_type, tarch_type, back_cast, T, 1)



%not much improvement with additional lags either p or q 


%Log likelihood on SP500 log returns , distribution innovation


% 17684  density garch

% 18447.4 normal with abs.value

% 18478 density scaled garch exp(-1/(2*(data(t-1)-theta(3))^2))

% 18477 normal

% 18488 engle lee (1,1)

% 18520 theta(3) * exp(-1/(2*(data(t-1)-0.19)^2))

% 18543 garch HN

% 18594 normal + asymmetric param 
% 18594.1 student t
% 18600.9 ged
% 18628.9 skew t
% 18749.8 skew t + asymmetric param



