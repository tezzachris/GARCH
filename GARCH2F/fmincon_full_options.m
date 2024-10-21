

%Numerical contrained optimization for the joint-likelihood of GARCH2F 

%Inputs:
%[ret]:[double Tx1] log returns
%[data]:[double Nx1] options dataframe - check loglik_option_full_iv to see what is needed
%[rf]:[double 1x1] risk free rate
%[prev]:[double px1] initial paramters, p=number of paramters

%Outputs:
%[a]:[double px1] MLEs
%[log_lik]:[double 1x1] maximum log-likelihood
%[exitflag]:[double 1x1] fmincon - optimizer exit flag
%[fo]:[double 1x1] fmincon - first order condition 

function [a,log_lik,start] = fmincon_full_options(ret, data, rf, prev)
    if isempty(prev)
       prev = start_full(ret,rf);
    end
    options = optimoptions(@fmincon, ...
                            'Display','off', ...
                            'Diagnostics','off', ...
                            'ConstraintTolerance',1e-7, ...
                            'algorithm','sqp');  %altri metodi interior point
    Aineq = []; 
    bineq = [];
    Aeq=[];
    beq=[];
    %lowerbound
    lb = [ 1e-50,1e-50,1e-50,-1000,1e-50,    1e-50,-1000,1e-50, 1e-50,1e-50,1e-50,1e-50,-100 ];
    %upperbound
    ub = [ 0.1,0.1,0.1,+1000,1,    0.1,+1000,1, 0.1,0.1,1,1,     +100 ];
    
    nonlcon = @(param) nonlincon_full(param);    
    estim_flag = 0; %flag for fminunc, in this case we use fmincon
    
    T = size(ret,1);
    N = size(data,1);
        
    fun=@(theta) (T+N)/(2*T) * loglik_full(theta,ret,estim_flag,rf) ... 
                   + (T+N)/(2*N) * loglik_option_full_iv(theta,data,estim_flag,rf) ;
    start = prev;
    [a,fval]=fmincon(fun,prev,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
    log_lik = -fval;

end 




