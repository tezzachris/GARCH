

%uncontrained optimization the likelihood to obtain the coefficients of GARCH
%model

function [a] = fmincon_hn(data)
    options = optimset('fmincon');
    options.Display = 'iter'; 
    options.TolX =  1e-20; 
    options.Tolfun = 1e-18;
    options.Maxiter = 1000; 
    options.MaxFunEvals =  5000; 
    options.HessUpdate = 'steep-dsc';  %use for problems with many variables 
    
    theta0 = [5.02*10^(-6) 1.32*10^(-6)  190.39  0.589 0.205];
    theta0 = [2.7184e-12     5e-6   150    0.8 3 ];%omega,alpha,gamma,beta,lambda1,h0
    A = []; 
    b = [];
    Aeq=[];
    beq=[];
    lb = zeros(1,length(theta0));
    ub = [];
    nonlcon = @nonlinearconstraints;
    fun=@(a)log_likelihood(data,a);
    a=fmincon(fun,theta0,A,b,Aeq,beq,lb,ub,nonlcon,options);
  
     %out=table(a(1),a(2),a(3),a(4),VariableNames={'omega';'beta';'alpha';'gamma'});
end 

function [c,ceq] = nonlinearconstraints(a)
    c = a(4) + a(2)*a(3)^2  - 1; %beta1+alpha1*gamma1^2 < 1
    ceq = [];
end

function val = log_likelihood(data, a) 

h = zeros(length(data),1); v=h; mu=h; z = h;
h(1) = ( a(1) + a(2) ) / ( 1-a(2)*a(3)^2-a(4) );
rf = 0;
mu(1) = rf + a(5)*h(1);
v(1) = - 0.5*( log(2*pi) +  log(h(1)) +  (data(1)-mu(1))^2 / h(1) );
z(1) = (data(1)-mu(1))/sqrt(h(1));

for t = 2:length(data)
    h(t) = a(1) + a(2) * ( z(t-1) - a(3)*sqrt(h(t-1))  )^2 + a(4) * h(t-1);
    mu(t) = rf + a(5)*h(t);
    z(t) = (data(t)-mu(t))/sqrt(h(t));
    v(t) = - 0.5 * ( log(2*pi) + log(h(t)) + (data(t)-mu(t))^2 / h(t) );
end
val = - sum(v); %we want to max so we put minus
end 

% function [r] = simul_hn (params,T)
%     eps=zeros(T,1); h=zeros(T,1);  
%     zeta= randn(0, 1, [T, 1]); %simulate zeta normal rv
%     omega = params(1);beta = params(2);alpha = params(3); gamma=params(4);
%     h(1) = ( omega+alpha ) / (1-alpha*gamma*gamma-beta); %initial volatility
%     eps(1)=sqrt(h(1))*zeta(1);
%         for i = 2:T
%             h(i) = omega + beta * h(i-1) + alpha * ( zeta(i-1)  - gamma*sqrt(h(i-1)) )^2  ;
%             eps(i) = sqrt(h(i))*zeta(i);
%             r
%         end
% end 
