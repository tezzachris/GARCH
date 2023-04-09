
%constrained likelihood optimization
function [a] = fmincon_garch(ret)

    options  =  optimset('fmincon');
    options  =  optimset(options , 'TolFun'      , 1e-20);
    options  =  optimset(options , 'TolX'        , 1e-10);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'off');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'Maxiter'      , 50000);
    options  =  optimset(options , 'MaxFunEvals' , 400*6);
    options = optimset(options, 'HessUpdate', 'bfgs');
    options = optimset(options, 'DiffMinChange', 0.00001);
    options = optimset(options, 'DiffMaxChange', 0.001);
    options = optimset(options, 'UseParallel', false);
  
    %theta0 = starting_values(data,1,0,1);
    p0 = [0.00001, 0.001, 0.8, 4];
    estim_flag=0;
    rf = 1e-7;
    p0=fminsearch('log_likelihood' , p0 , optimset('display','iter','Maxiter',1000),ret,estim_flag,rf);

    A = []; %beta + alpha < 1
    b = [];
    Aeq=[];
    beq=[];
    lb = [0, 0, 0, 0];
    ub= [];
    fun = @(p0) log_likelihood(p0,ret,estim_flag,rf);
    nonlcon = 'nonlinearconstraints';

    a=fmincon(fun,p0, A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    %se=sqrt(diag(inv(hessian))); %standard errors
    
    %MODEL COMPARISON
    %aic = 2 * length(a) - 2*log_lik; 
    %bic= - 2 * log_lik + length(a)*log(length(data)); %higher is worse
    
    %rmse = sqrt (  sum( (data-fit).^2 ) / length(data) ) ;

end

