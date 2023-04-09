%Christoffersen 2008 
function [a,log_lik , exitflag, fo] = fminunc_chris2008(ret,rf,prevparam)

    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-10);
    options  =  optimset(options , 'TolX'        , 1e-20);
    options  =  optimset(options , 'Display'     , 'off');
    options  =  optimset(options , 'Diagnostics' , 'off');
    options  =  optimset(options , 'LargeScale'  , 'off');
    nparams = 8;
    options  =  optimset(options , 'MaxFunEvals' , 200*nparams);

    estim_flag = 1;
    if isempty(prevparam)
        prevparam = startingvals_chris(ret,rf);
    end

    %Transform to Unconstrained parameters
    theta = paramtotheta(prevparam); 

    % parameter estimation using the custom function log_lik & the MFE version
    [theta,fval,exitflag,output] = fminunc('log_likelihood_chris', theta , options, ret, estim_flag,rf);
    
% transforming back the parameters since the optimization returns theta
    a = thetatoparam(theta);
    log_lik = -fval;
    fo=output.firstorderopt;
    [VCV] = hessian_2sided(a,ret,rf);
     se = sqrt(diag((VCV))); %standard errors
     signif = ( abs( a ) ./ se' ) >= tinv(0.975, length(ret)-length(a)) ;
     aic = 2 * length(a) - 2*log_lik; 
     bic= - 2 * log_lik + length(a)*log(length(ret)); %higher is worse

  
end 