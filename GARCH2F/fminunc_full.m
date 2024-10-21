function [a, log_lik , exitflag, fo] = fminunc_full(ret,rf ,prevparam)

    options  =  optimset('fminunc');
    options  =  optimset(options , 'Display'     , 'off');
    options  =  optimset(options , 'Maxiter'      , 50000);
    nparams = 13;
    options  =  optimset(options , 'MaxFunEvals' , 400*nparams);

    estim_flag = 1;
    
    if isempty(prevparam)
        prevparam = start_full(ret,rf);
    end

    theta = paramtotheta_full(prevparam);

    fun=@(theta) loglik_full(theta,ret,estim_flag,rf) ;
    
    [theta,fval,exitflag,output] = fminunc( fun , theta, options );
    
    a = thetatoparam_full(theta);
    log_lik = -fval;
    fo=output.firstorderopt;

end 
