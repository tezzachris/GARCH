
%Estimate parameters for GARCH(1,1) in the following order
%omega, alpha, beta, lambda
function [a, log_lik , exitflag, fo] = fminunc_garch(ret,rf,prevparam)

    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-20);
    options  =  optimset(options , 'TolX'        , 1e-10);
    options  =  optimset(options , 'Display'     , 'off');
    options  =  optimset(options , 'Diagnostics' , 'off');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'Maxiter'      , 50000);
    options  =  optimset(options , 'MaxFunEvals' , 400*6);
    options = optimset(options, 'HessUpdate', 'bfgs');
    options = optimset(options, 'DiffMinChange', 0.00001);
    options = optimset(options, 'DiffMaxChange', 0.001);
    options = optimset(options, 'UseParallel', false);

    estim_flag = 1;
    if isempty(prevparam)
        theta = startingvals(ret,rf);
    else
        theta = paramtotheta(prevparam);
    end
    %param is Constrained params  
    %theta is Unconstrained parameters

    [theta,fval,exitflag,output] = fminunc( 'log_likelihood', theta , options, ret, estim_flag, rf);
    
    %transforming back the parameters since the optimization returns theta
    a = thetatoparam(theta);
    log_lik = -fval;
    fo=output.firstorderopt;

     [VCV] = hessian_2sided(a,ret,rf);
     se = sqrt(diag((VCV))); %standard errors
     signif = ( abs( a ) ./ se' ) >= tinv(0.975, length(ret)-length(a)) ;
     aic = 2 * length(a) - 2*log_lik; 
     bic= - 2 * log_lik + length(a)*log(length(ret)); %higher is worse
   
    
end 