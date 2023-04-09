function [ rt ] = vecsimulate_hn(T,a,rf,nsim)
    
    h1=zeros(T,nsim);  mu=zeros(T,nsim);  
    z1 = randn(T, nsim); 
    h1(1)=(a(1)+a(2)) / (1 - a(2)*a(3)^2 - a(4));
    e1(1,:)= sqrt(h1(1,:)).*z1(1,:); 
    mu(1,:)= rf + a(5) * h1(1,:);

    for t = 2:T
        h1(t,:) = a(1) + a(2) * (z1(t-1,:)  - a(3).*sqrt(h1(t-1,:))).^2  + a(4) * h1(t-1,:);
        e1(t,:) = sqrt(h1(t,:)).*z1(t,:); 
        mu(t,:) = rf + a(5) * h1(t,:); 
    end
     rt = mu + e1;
end 