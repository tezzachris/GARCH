% function [val,v] = log_likelihood(theta,data,estim_flag,rf) %theta are the transformed params
% 
%     if estim_flag == 1
%         a = thetatoparam(theta);
%     elseif estim_flag == 0
%         a = theta;
%     end
%    
%     q = zeros(length(data), 1); s = q;
%     q(1) = a(1)/(1-a(4)) ;
%     s(1) = 0 ;
%     v=zeros(length(data), 1); mu = v;
%     mu(1) = rf + a(8) * (q(1)+s(1));
%     v(1)= - 0.5 * ( log(2*pi) + log(q(1)+s(1)) + (data(1)-mu(1))^2 / (q(1)+s(1)) );
%     
%     for t = 2:length(data)
%         zeta = ( data(t-1)-mu(t-1) )/sqrt(q(t-1)+s(t-1));
%         q(t) = a(1) +  a(2) * ( zeta^2 - 1 - 2*a(3)*sqrt(s(t-1)+q(t-1))*zeta ) + a(4) * q(t-1);
%         s(t) = a(5) * ( zeta^2  - 1 - 2*a(6)*sqrt(s(t-1)+q(t-1))*zeta ) + a(7) * s(t-1);
%         mu(t) = rf + a(8) * (q(t)+s(t));
%         v(t) = - 0.5 * ( log(2*pi) + log(q(t)+s(t)) + (data(t)-mu(t))^2 /(q(t)+s(t)) );
%     end
%     val = - sum(v); %we want to max so we put minus
%     v = - v ; 
% end 

function [val,v] = log_likelihood_chris(theta,data,estim_flag,rf) %theta are the transformed params

    if estim_flag == 1
        a = thetatoparam(theta);
    elseif estim_flag == 0
        a = theta;
    end
   
    qt = zeros(length(data), 1); ht = qt;
    qt(1) = a(1)/(1-a(4)) ;
    ht(1) = qt(1) ;
    v=zeros(length(data), 1); mu = v;
    mu(1) = rf + a(8) * ht(1);
    v(1)= - 0.5 * ( log(2*pi) + log(ht(1)) + (data(1)-mu(1))^2 / (ht(1)) );
    
    for t = 2:length(data)
        zt = ( data(t-1)-mu(t-1) )/sqrt(ht(t-1));
        qt(t) = a(1) +  a(2) * ( zt^2 - 1 - 2*a(3)*sqrt(ht(t-1))*zt ) + a(4) * qt(t-1);
        ht(t) = qt(t) + a(5) * ( zt^2  - 1 - 2*a(6)*sqrt(ht(t-1))*zt ) + (a(5)*a(6)^2 + a(7)) * (ht(t-1)-qt(t-1));
        mu(t) = rf + a(8) * (ht(t));
        v(t) = - 0.5 * ( log(2*pi) + log(ht(t)) + (data(t)-mu(t))^2 /(ht(t)) );
    end
    val = - sum(v); %we want to max so we put minus
    v = - v ; 
end 

