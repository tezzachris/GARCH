

%use risk neutral params for pricing
function [ logS, h1, rt ] = simulate(T,a,rf,St0)
    h1=zeros(T,1);  mu=h1;  
    z1 = randn(T, 1); 
    h1(1)= (a(1) + a(2))/(1-a(2)*a(3)^2-a(4)); %mean stationary process
    mu(1)= rf - a(5) * h1(1);
    logS = zeros(T,1); 
    logS(1)=log(St0);
    for t = 2:T
        h1(t) = a(1) + a(2) * ( z1(t-1) - a(3)*sqrt(h1(t-1)) )^2  + a(4) * h1(t-1);
        mu(t) = rf + a(5) * h1(t); 
        logS(t) = logS(t-1) + mu(t) + sqrt(h1(t))*z1(t);
    end
    rt = mu + sqrt(h1) .* z1  ;
end 