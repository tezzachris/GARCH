
function [res] = christoff1d(data)
    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-005);
    options  =  optimset(options , 'TolX'        , 1e-005);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');

% transformations of the initial values by using the inverse of the functions
    %p=1; o = 1; q=1;
    %params = starting_values(data, p, o, q ); %starting omega beta alpha
    %params = [0.032 ; 0.001 ; 0.9; 0.8 ]; %omega beta alpha
    
    p = array2table(  ...
        [0.118	, ...
        0.005 , 11.98   ,0.08, ...
        0.002 ,5, 0.06],...
               'VariableNames',{'omega'; 'alpha1';'gamma1';'beta1';'alpha2';'gamma2';'beta2'}); 
 
    options  =  optimset(options , 'MaxFunEvals' , 200*(2+size(p,2)));

    %OMEGA
    tomega = log(p.omega); 
    talpha1 = log(p.alpha1 ./ (1-p.alpha1));
    tgamma1 = log(p.gamma1); 
    tbeta1 = log(p.beta1./(1-p.beta1));
    talpha2 = log(p.alpha2 ./ (1-p.alpha2));
    tgamma2 = log(p.gamma2); 
    tbeta2 = log(p.beta2./(1-p.beta2));
    tparams = [tomega; talpha1;tgamma1; tbeta1;talpha2;tgamma2;tbeta2];
    
    fminsearch(@loglik_christoff , tparams , optimset('display','iter','Maxiter',2000),data);

% parameter estimation using the custom function log_lik & the MFE version
    [theta, fval, exitflag, output, grad, hessian] = fminunc( @loglik_christoff, tparams, options,data);
    
% transforming back the parameters since the optimization returns theta
    
    a(1)= exp(theta(1)); 
    a(2)=(exp(theta(2))/(1+exp(theta(2))));
    a(3)=exp(theta(3));
    a(4)=(exp(theta(4))/(1+exp(theta(4))));
    a(5)=(exp(theta(5))/(1+exp(theta(5))));
    a(6)=exp(theta(6));
    a(7)=(exp(theta(7))/(1+exp(theta(7))));
    
    log_lik = -fval; %log likelihood value at convergence 
    first_order  = output.firstorderopt; %first order optimality measure
    se=sqrt(diag(inv(hessian))); %standard errors
    n=length(data);
    res = [a;table2array(p)];
%     res=array2table([a a(4)+a(2)*a(3)^2 ; 
%                      se' , nan ; 
%                      transpose( abs( a' ) ./ se  >= tinv(0.975, n-length(a)) ), nan; %df = n - p
%                      params' , nan],...
%         'VariableNames',{'omega';'alpha';'gamma';'beta'; 'alpha*gamma^2+beta'},...
%         'RowNames',{'Heston-Nandi(1,1)';'Std.err.';'Signif 95%';'Initial_params'});
%     

end 

function val = loglik_christoff(theta,data) %data is zeta N(0,1) vector Tx1
              
    a(1)= exp(theta(1)); 
    a(2)=(exp(theta(2))/(1+exp(theta(2))));
    a(3)=exp(theta(3));
    a(4)=(exp(theta(4))/(1+exp(theta(4))));
    a(5)=(exp(theta(5))/(1+exp(theta(5))));
    a(6)=exp(theta(6));
    a(7)=(exp(theta(7))/(1+exp(theta(7))));
    
    h1 = zeros(length(data), 1); h2= h1;
    h1(1) = a(1)/(1-a(4)); 
    h2(1) = 0;

    v=zeros(length(data), 1);
    v(1)= - 0.5* ( log(2*pi) + log(h1(1)+h2(1)) + (data(1))^2 /(h1(1)+h2(1)) );
    
    for t = 2:length(data)

        zeta = data(t-1) / sqrt(h1(t-1)+h2(t-1));

        h1(t) = a(1) + a(2) * ( zeta^2 - 1 - 2*a(3)*sqrt(h1(t-1)+h2(t-1))*zeta ) + a(4) * h1(t-1) ;
        h2(t) = a(7)  * h2(t-1) + a(5) * ( zeta^2 - 1 - 2*a(6)*sqrt(h1(t-1)+h2(t-1))*zeta  );
            
        v(t) = -0.5 * (log(2*pi) + log(h1(t)+h2(t)) + data(t)^2 / (h1(t)+h2(t)) );
    end
val = - sum(v); %we want to max so we put minus
end 


function [eps,seed] = simul_garch
    seed=floor(rand(1)*100000);
    rng(seed);
    T=5000; 
    eps=zeros(T,1); h=zeros(T,1);  
    zeta= randn(T,1); 
    omega = 0.003; 
    alpha = 0.01;
    beta = 0.96;
    h(1) = (omega)/(1-beta-alpha); %initial volatility
    eps(1)=sqrt(h(1))*zeta(1);
        for i = 2:T
            h(i) = omega + beta * h(i-1) + alpha * eps(i-1).^2;
            eps(i) = sqrt(h(i))*zeta(i);
        end
end 



function [startingvals]=starting_values(data, p, o, q )

% Perform a grid search to find decent starting values for TARCH(P,O,Q)
% esimtation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on error_type.
%
% USAGE:
%   [STARTINGVALS,NU,LAMBDA,LLS,OUTPUT_PARAMETERS] = ...
%        tarch_starting_values(STARTINGVALS,DATA,FDATA,FIDATA,P,O,Q,T,ERROR_TYPE,TARCH_TYPE);
%
% INPUTS:
%   STARTINGVALS     - A vector of starting values or empty to perform a grid search
%   DATA             - Vector of mean zero residuals
%   FDATA            - Either abs(data) or data.^2, depending on tarch_type
%   FIDATA           - fdata times an indicator for negative, e.g. fdata.*(data<0)
%   P                - The lag order length for ARCH
%   O                - The lag order of asymmetric terms
%   Q                - The lag order length for GARCH
%   T                - Length of data
%   ERROR_TYPE       - The type of error being assumed, valid types are:
%                        1 if 'NORMAL'
%                        2 if 'STUDENTST'
%                        3 if 'GED'
%                        4 if 'SKEWT'
%   TARCH_TYPE        - 1 for absolute vale of return
%                     - 2 for squared returns (standard case)
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (1+p+o+q) by 1
%   NU                - Distribution kurtosis parameter, empty if not applicable
%   LAMBDA            - Distribution asymmetry parameter, empty if not applicable
%   LLS               - A vector of log likelihoods corresponding to OUTPUT_PARAMETERS
%   OUTPUT_PARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%Initialize variables
    LLs=[];
    output_parameters=[];


    %Procedure is to find best starting values, using a grid search find values for normal, then

    %Possible starting values based on commonly estimated values
    a=[.05 .1 .2]; %omega
    la=length(a);

    g=[.01 .05 .2]; %alpha
    lg=length(g)+1;

    agb=[.5 .8 .9 .95 .99]; %beta
    lb=length(agb); 
    
    %Many output_parameters and LLs
    output_parameters=zeros(la*lb*lg,1+p+o+q);
    LLs=zeros(la*lb*lg,1);

    %Adjustment is needed to the intercept.  Assumes normality
    adj_factor=1;
    back_cast=cov(data);
    
    covar=cov(data);

    %Use an index to count
    index=1;

    for i=1:la
        %Loop over omega 
        alpha=a(i);
        for j=1:lg
            %Loop over alpha
            if j==lg
                gamma=-alpha/2;
            else
                gamma=g(j);
            end

            for k=1:lb
                %Loop over beta

                temp_alpha=alpha;
                temp_gamma = [];
                %Make sure gamma satisfies necessary constraints
                %Beta must also satisfy the same constraints
                beta=agb(k)-sum(temp_alpha)-0.5*sum(temp_gamma);
                %Pick omega to match the unconditional
                omega=covar*(1-sum(temp_alpha)*adj_factor-0.5*sum(temp_gamma)-sum(beta));
     
                %Build the parameter vector
                parameters= [omega; temp_alpha; beta];
             
                %Set the output parameters
                output_parameters(index,:)=parameters';
                parameters=[parameters(1);parameters(3);parameters(2)];
                %Set the log likelihoods
                LLs(index) = log_lik_garch(parameters, data);
                %Increment
                index=index+1;
            end
        end
    end
    %Sort the LLs so the best (lowest, since we minimize the -1*LL)
    [LLs,index]=sort(LLs);
    %Order the starting values
    startingvals=output_parameters(index(1),:)';
    %Order the ouput_parameters
    output_parameters=output_parameters(index,:);

end 



%P.474 long run and short run component vol model engle lee

%rt  = mut + et, mut= rf + lambda1 * qt + lambda2 * st
%et  = sqrt(ht) * zt, zt is N(0,1)
%ht  = qt + st
%qt = omega  + rho * qt-1 + alfa1 * ( et-1 ^ 2 - ht-1 )
%st  = beta * st-1 + alpha2 * ( et-1 ^ 2 - ht-1 )

% function [ params , compare  ] = mle_christoff(data)
%    
%    options = optimset('TolX', 1e-20, 'TolFun', 1e-10, 'Maxiter', 50000, 'MaxFunEvals', 50000,...
%                         'HessUpdate', 'bfgs', 'UseParallel', true, 'display', 'iter');
%     %rng('default') %for reproducibility
%     %T=5000; 
%     theta0=[ 8e-4; 0.989; 2e-4; 6e1; %omega; rho; alpha1; gamma1;
%                    0.7; 3e-5; 5e1; % beta; alpha2; gamma2;
%                     1 ];  %lambda1; 
% 
%     A = []; % No other constraints
%     b = [];
%     Aeq = [];
%     beq = [];
%     lb = [0; 0; 0; -inf; 0; 0; -inf; -inf ]; % contraints for positive cond variance
%     ub = [];
%     nonlcon=@nonlinearconstraints;
%     fun=@(theta0)log_lik_2factor_garch(data,theta0);
%     
%     [a,fval,exitflag,output,lambda,grad,hessian]  = fmincon(fun, theta0, A, b ,Aeq,beq,lb,ub, nonlcon ,options);
%     
%     se=sqrt(diag(inv(hessian))); %standard errors
% 
%     params = array2table([a'; ...
%                     se';...
%                     transpose( abs( a ./ se ) >= tinv(0.95,length(data)-length(a)) ); ...
%                     theta0'],...
%         'VariableNames',{'omega';'rho';'alpha1';'gamma1';'beta';'alpha2';'gamma2';'lambda1'},...
%         'RowNames',{'Christoff';'Std.err.';'Signif. 90%';'Initial'});
%     
%     aic= 2*length(a) - 2*fval; bic= - 2*fval + length(a)*log(length(data));
%     
%     [fit]=simul_2factor_garch(a,length(data));
%     rmse = sqrt (  sum( (data-fit).^2 ) / length(data) ) ;
%     compare = array2table([ fval, aic, bic , rmse ,length(data) ], 'VariableNames',{'Log-lik';'AIC';'BIC';'RMSE';'Observations'});
%     
%     
%   end 
% 
% function [c,ceq] = nonlinearconstraints(a)
%     c(1) = a(2)  - 1; %rho <= 1
%     c(2) = a(5) + a(6)*a(7)^2 - a(2)  ; %beta+alpha2 <= rho
%     ceq = [];
% end
% 
% function [ val ] = log_lik_2factor_garch(data,a) 
%     
%     h1 = zeros(length(data), 1); h1(1) = a(1)/(1-a(2)); 
%     h2 = zeros(length(data), 1); h2(1) = 0;
%     rf=0.02/100;
%     mu = zeros(length(data), 1); mu(1) = rf + a(8) * (h1(1)+ h2(1));
%     v = [zeros(length(data), 1)]; 
%     v(1) = 1/2 * ( log(2*pi) + log(h1(1)+h2(1)) + ( data(1) - mu(1) )^2 / ( h1(1)+h2(1) ));
%     
%     for t = 2:length(data)
%             z1(t-1) = ( data(t-1) - mu(t-1) ) / sqrt(h1(t-1)+h2(t-1));
%             h1(t) = a(1) + a(2) * h1(t-1) + a(3) * ( z1(t-1)^2 - 1 - 2*a(4)*sqrt(h1(t-1)+h2(t-1))*z1(t-1)  ) ;
%             h2(t) = (a(5) + a(6)*a(7)^2 ) * h2(t-1) + a(6) * ( z1(t-1)^2 - 1 - 2*a(7)*sqrt(h1(t-1)+h2(t-1))*z1(t-1)  );
%             mu(t) = rf + a(8) * (h1(t) +  h2(t));
%             v(t) = log(2*pi) + log(h1(t)+h2(t))  + ( data(t) - mu(t) )^2 / ( h1(t)+h2(t) ) ;
%     end
%     val = 1/2*sum(v); 
% end 
% 
%     
% function [ r ] = simul_2factor_garch(a,T)
%     e1=zeros(T,1); h1=zeros(T,1); h2=zeros(T,1);    
%     z1= randn(T, 1); %simulate zeta normal rv
%     h1(1)=a(1)/(1-a(2)); %mean stationary process
%     h2(1)=0; 
%     e1(1) = sqrt(h1(1)+h2(1))*z1(1); 
%     r = zeros(T,1); mu = zeros(T,1); 
%     rf = 0.02/100;
%     mu(1) = rf + a(8) * (h1(1) + h2(1));
%     r(1) = mu(1) + e1(1); 
%         for t = 2:T
%             %Long process qt
%             h1(t) = a(1) + a(2) * h1(t-1) + a(3) * ( z1(t-1)^2 - 1 - 2*a(4)*sqrt(h1(t-1)+h2(t-1))*z1(t-1)  ) ;
%             %Short process st
%             h2(t) = (a(5) + a(6)*a(7)^2 ) * h2(t-1) + a(6) * ( z1(t-1)^2 - 1 - 2*a(7)*sqrt(h1(t-1)+h2(t-1))*z1(t-1)  );
%             %Innovation
%             e1(t) = sqrt(h1(t)+h2(t))*z1(t); 
%             mu(t) = rf + a(8) * (h1(t) + h2(t));
%             r(t) = mu(t) + e1(t);
%         end
% 
% end 


