function [ rt ] = simulate(T,a,rf)
    
    e1 = zeros(T,1); h1=zeros(T,1);  mu=zeros(T,1);  
    z1 = randn(T, 1); 
    h1(1)= a(1)/(1-a(2)-a(3)); %mean stationary process
    e1(1)= sqrt(h1(1))*z1(1); 
    mu(1)= rf + a(4) * h1(1);

    for t = 2:T
        h1(t) = a(1) + a(2) * e1(t-1)^2  + a(3) * h1(t-1);
        e1(t) = sqrt(h1(t))*z1(t); 
        mu(t) = rf + a(4) * h1(t); 
    end
     rt = mu + e1 ;
end 