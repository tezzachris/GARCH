
%P.474 long run and short run component vol model engle lee

%rt  = mut + et, mut= rf + lambda1 * qt + lambda2 * st
%et  = sqrt(ht) * zt, zt is N(0,1)
%ht  = qt + st
%qt = omega  + rho * qt-1 + alfa1 * ( et-1 ^ 2 - ht-1 )
%st  = beta * st-1 + alpha2 * ( et-1 ^ 2 - ht-1 )

function [ params , compare  ] = mle_engle_lee(data)
   
   options = optimset('TolX', 1e-20, 'TolFun', 1e-5, 'Maxiter', 50000, 'MaxFunEvals', 50000,...
                        'HessUpdate', 'bfgs', 'UseParallel', true, 'display','iter');
   
    theta0=starting_values(data);  %lambda1; lambda2

    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [0; 0; 0; 0; 0]; % contraints for positive cond variance
    ub = [];
    nonlcon=@nonlinearconstraints;
    fun=@(theta0)log_lik_2factor_garch(data,theta0);
    
    [a,fval,exitflag,output,lambda,grad,hessian]  = fmincon(fun, theta0, A, b ,Aeq,beq,lb,ub, nonlcon ,options);
    
    log_lik = -fval;
    se=sqrt(diag(inv(hessian))); %standard errors

    params = array2table([a'; ...
                          se';...
                          transpose( abs( a ) ./ se  >= tinv(0.975,length(data)-length(a)) ); ...
                          theta0'],...
        'VariableNames',{'omega';'rho';'alpha1';'beta';'alpha2'},...
        'RowNames',{'Engle-Lee 1999 (1,1)';'Std.err.';'Signif. 95%';'Initial params'});
    
    aic= 2*length(a) - 2*log_lik; 
    bic= - 2*log_lik + length(a)*log(length(data));
    
    %rmse = sqrt (  sum( (data-fit).^2 ) / length(data) ) ;
    compare = array2table([ log_lik, aic, bic ,length(data) ], 'VariableNames',{'Log-lik';'AIC';'BIC';'Observations'});
    
    
  end 

function [c,ceq] = nonlinearconstraints(a)
    c(1) = a(2)  - 1; %rho <= 1
    c(2) = a(4) + a(5)  - a(2)  ; %beta+alpha2 <= rho
    %c(3) = - a(4) - a(5); % beta + alpha2 > 0
    %c(4) = a(3) - a(5); %alpha1 < alpha2 (ie st more reactive to shocks than qt)
    ceq = [];
end

function [ val ] = log_lik_2factor_garch(data,a) 
    
    h1 = zeros(length(data), 1); h1(1) = a(1)/(1-a(2)); 
    h2 = zeros(length(data), 1); h2(1) = 0;
    rf=0.02/100;
    mu = zeros(length(data), 1); mu(1) = 0; % rf + a(6) * h1(1) + a(7) * h2(1);
    v = [zeros(length(data), 1)]; 
    v(1) = - 1/2 * ( log(2*pi) + log(h1(1)+h2(1)) + ( data(1) - mu(1) )^2 / ( h1(1)+h2(1) ) );
    
    for t = 2:length(data)
            h1(t) = a(1) + a(2) * h1(t-1) + a(3) * ( (data(t-1)-mu(t-1) )^2 - (h1(t-1)+h2(t-1))  ) ;
            h2(t) = a(4) * h2(t-1) + a(5) * ( (data(t-1)-mu(t-1) )^2 - (h1(t-1)+h2(t-1)) ) ;
            mu(t) = 0; %,rf + a(6) * h1(t) + a(7) * h2(t);
            v(t) = - 0.5 * ( log(2*pi) + log(h1(t)+h2(t))  + ( data(t) - mu(t) )^2 / ( h1(t)+h2(t) ) );
    end
    val = - sum(v); 
end 

    
function [ r ] = simul_2factor_garch(a,T)
    e1=zeros(T,1); h1=zeros(T,1); h2=zeros(T,1);    
    z1= randn(T, 1); %simulate zeta normal rv
    h1(1)=a(1)/(1-a(2)); %mean stationary process
    h2(1)=0; 
    e1(1) = sqrt(h1(1)+h2(1))*z1(1); 
    r = zeros(T,1); mu = zeros(T,1); 
    rf = 0.02/100;
    mu(1) = 0; %rf + a(6) * h1(1) + a(7) * h2(1);
    r(1) = mu(1) + e1(1); 
        for t = 2:T
            %Long process qt
            h1(t) = a(1) + a(2) * h1(t-1) + a(3) * ( e1(t-1)^2 - (h1(t-1)+h2(t-1))  ) ;
            %Short process st
            h2(t) = a(4) * h2(t-1) + a(5) * ( e1(t-1)^2 - (h1(t-1)+h2(t-1)) ) ;
            %Innovation
            e1(t) = sqrt(h1(t)+h2(t))*z1(t); 
            mu(t) = 0; %rf + a(6) * h1(t) + a(7) * h2(t);
            r(t) = mu(t) + e1(t);
        end

end 

function [startingvals]=starting_values(data)

% Perform a grid search to find decent starting values for TARCH(P,O,Q)
% esimtation.  If starting values are user supplied (and thus nonempty), reformats
% starting values depending on error_type.

%Initialize variables
   

    %Procedure is to find best starting values, using a grid search find values for normal, then

    %Possible starting values based on commonly estimated values
   
    alpha1=[.01 .05 .1]; %alpha
    lalpha1=length(alpha1)+1;

    rho=[.5 .8 .9 .95 .99]; %beta
    rho = rho' - alpha1;
    rho = reshape(rho, size(rho,1)*size(rho,2),1) ; %convert into vector
    lrho=length(rho); 
    
    %omega=[.05 .1 .2]; %omega
    omega = cov(data) * ( 1 - rho );
    
    log_lik = zeros( length(omega) , 1, length(alpha1));
    

    for i = 1:length(alpha1)
        %mega_rho = nchoosek( [omega ; rho] , 2);
        parameters = [omega ,  rho, ...
                      alpha1(i)*ones(length(omega),1), ... 
                      rho - 5/1000, ...
                      alpha1(i)*ones(length(omega),1) + 5/10000]; %alpha2  
        log_lik(:,:,i) =  arrayfun(@(ROWIDX) log_lik_2factor_garch(data,parameters(ROWIDX,:)), (1:size(parameters,1)).');
    end
    
    log_lik = cat(3,log_lik);
    [log_lik,index]=sort(log_lik);
    startingvals = parameters(index(1),:)';

end
