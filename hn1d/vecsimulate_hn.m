function [ rt, h1, z1 ] = vecsimulate_hn(T,a,rf,nsim)
    
    h1 = zeros(T,nsim);  
    z1 = randn(T, nsim); 
  
    for t = 1:T
        if t == 1
            h1(t,:)= (a(1)+a(2)) / (1 - a(2)*a(3)^2 - a(4));
        else
            h1(t,:) = a(1) + a(2) * ( z1(t-1,:)  - a(3).*sqrt(h1(t-1,:))).^2  + a(4) * h1(t-1,:);        
        end    
    end
    rt = rf + a(5) * h1 + sqrt(h1).*z1;
end 