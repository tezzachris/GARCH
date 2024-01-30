

%Numerical Contrained optimization for the likelihood of GARCH(1,1) Heston
%Nandi type

%Inputs:
%ret:[double nx1] vector of log returns
%rf:[double 1x1] risk free rate

%Outputs:
%Coefficients that maximize the likelihood

function [a,log_lik,exitflag,fo] = fmincon_hn(ret, rf, prev)
    if isempty(prev)
       prev = startingvals_hn(ret,rf);
    end
    options = optimoptions(@fmincon,'Display','off','MaxFunEvals',5*10000); %altri metodi interior point
    Aineq = []; 
    bineq = [];
    Aeq=[];
    beq=[];
    lb = [1e-15    1e-15   -1000  0.1 -100];
    ub = [0.1 0.1   1000  1   +100];
    nonlcon = @(param) nonlinearconstraints(param);    
    estim_flag = 0;
    fun=@(theta) log_likelihood_hn(theta,ret,estim_flag,rf);
    [a,fval,exitflag,fo]=fmincon(fun,prev,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
    log_lik = -fval;

end 




