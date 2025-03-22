



function [a,log_lik,flag,fo] = fmincon_full(ret, rf, prev)
    
    if isempty(prev)
        prev = start_full(ret,rf);
    end
    
    options = optimoptions(@fmincon,'Display','none');
    %lowerbound
    lb = [-1e-12,-1e-12,-1e-12,-1000,-1e-12,    -1e-12,-1000,-1e-12,    -1e-12,-1e-12, -1e-12,-1e-12,   -5];
    %upperbound
    ub = [0.1,0.1,0.1,+1000,1,    0.1,+1000,1,   0.1,0.1,  1,1,   +5];
    Aineq = []; 
    bineq = [];
    Aeq=[];
    beq=[];
    nonlcon = @(param) nonlincon_full(param);  
    estim_flag = 0;
    fun=@(param) loglik_full(param,ret,estim_flag,rf);
    [a,fval,flag,fo]=fmincon(fun,prev,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
    log_lik = -fval;
end 




