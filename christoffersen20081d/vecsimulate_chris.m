% function [ rt ] = vecsimulate(T,a,rf,nsim)
%     
%     h1=zeros(T,nsim); h2=h1;  mu=zeros(T,nsim);  
%     z1 = randn(T, nsim); 
%     h1(1,:)= a(1)/(1-a(4)); %mean stationary process
%     h2(1,:)= 0;
%     mu(1,:)= rf + a(8) * (h1(1,:)+h2(1,:));
% 
%     for t = 2:T
%         h1(t,:) = a(1) + a(2) * ( z1(t-1,:).^2 - 1 - 2*a(3)*sqrt(h1(t-1,:)+h2(t-1,:)).*z1(t-1,:) )  + a(4) * h1(t-1,:);
%         h2(t,:) = a(5) * ( z1(t-1,:).^2 - 1 - 2*a(6)*sqrt(h1(t-1,:)+h2(t-1,:)).*z1(t-1,:) )  + (a(5)*a(6)^2+a(7)) * h2(t-1,:); 
%         mu(t,:) = rf + a(8) * (h1(t,:)+h2(t,:)); 
%     end
%     rt = mu + z1.*sqrt(h1+h2) ;
% end 


function [ rt ] = vecsimulate_chris(T,a,rf,nsim)
    
    qt=zeros(T,nsim); ht=qt;  mu=qt;  
    zt = randn(T, nsim); 
    qt(1,:)= a(1)/(1-a(4)); %mean stationary process
    ht(1,:)= qt(1,:);
    mu(1,:)= rf + a(8) * ht(1,:);

    for t = 2:T
        qt(t,:) = a(1) + a(2) * ( zt(t-1,:).^2 - 1 - 2*a(3)*sqrt(ht(t-1,:)).*zt(t-1,:) )  + a(4) * qt(t-1,:);
        ht(t,:) = qt(t,:) + a(5) * ( zt(t-1,:).^2 - 1 - 2*a(6)*sqrt(ht(t-1,:)).*zt(t-1,:) )  + ( a(5)*a(6)^2 + a(7) ) * ( ht(t-1,:)-qt(t-1,:) ); 
        mu(t,:) = rf + a(8) * ht(t,:); 
    end
    rt = mu + zt.*sqrt(ht) ;
end 

