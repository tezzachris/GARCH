

%Numerical Contrained optimization for the likelihood of GARCH(1,1) Heston
%Nandi type

%Inputs:
%ret:[double nx1] vector of log returns
%rf:[double 1x1] risk free rate

%Outputs:
%Coefficients that maximize the likelihood

function [loglik] = fmincon_hn2d(ret, rf, filter)
    rng default % For reproducibility
    options = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); %altri metodi interior point

    ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',false);

    gs = GlobalSearch(ms);%gs = GlobalSearch;

    %theta0 = [5.02*10^(-6) 1.32*10^(-6)  190.39  0.589 0.205];
    %theta0 = [2.7184e-12     5e-6   150    0.8 3 ];%omega,alpha,gamma,beta,lambda1,h0
    param = startingvals_hn2d(ret,rf,filter);

    A = []; 
    b = [];
    Aeq=[];
    beq=[];
    lb = [0 0 0 0 0 0 0 -inf -inf]; %lambda can be negative
    ub = [];
    nonlcon = @(param) nonlinearconstraints(param); 
    estim_flag = 0;

    fun=@(param) log_likelihood_hn2d(param,ret,estim_flag,rf,filter);
    %a=fmincon(fun,param,A,b,Aeq,beq,lb,ub,nonlcon,options);

    problem = createOptimProblem('fmincon','x0',param,...
        'objective',fun,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
    
    [a,fval] = run(gs,problem);
    
    loglik = -fval; 
  
end 




