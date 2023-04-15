
function [val,v] = log_likelihood_hn2d(theta,data,estim_flag,rf,filter) 
  
    if estim_flag == 1
        a = thetatoparam(theta);
    elseif estim_flag == 0
        a = theta;
    end

h1 = zeros(length(data), 1); v=h1; mu = h1; e1=h1; e2=h1; h2 = h1;

for t = 1:length(data)
    if t==1
        h1(t)=(a(1)+a(2)) / (1 - a(2)*a(3)^2 - a(4));
        h2(t)=(a(5)/(1-a(7)-a(5)*a(6)^2));
    else 
        if filter == 1
            e1(t-1) = sqrt(h1(t-1)/(h1(t-1)+h2(t-1)))*(data(t-1)-mu(t-1));
            e2(t-1) = data(t-1) - mu(t-1) - e1(t-1);
        elseif filter == 2
            e2(t-1) = sqrt(h2(t-1)/(h1(t-1)+h2(t-1)))*(data(t-1)-mu(t-1));
            e1(t-1) = data(t-1) - mu(t-1) - e2(t-1);
        elseif isempty(filter)
            e1(t-1) = sqrt(h1(t-1)/(h1(t-1)+h2(t-1)))*(data(t-1)-mu(t-1));
            e2(t-1) = sqrt(h2(t-1)/(h1(t-1)+h2(t-1)))*(data(t-1)-mu(t-1));
        end
                    
        h1(t) = a(1) + a(2) * ( e1(t-1)/sqrt(h1(t-1)) - a(3)*sqrt(h1(t-1)) )^2 + a(4) * h1(t-1);   
        h2(t) = a(5) * ( e2(t-1)/sqrt(h2(t-1)) - a(6)*sqrt(h2(t-1)) )^2 + a(7) * h2(t-1);
    end
    mu(t) = rf + a(8)*h1(t) + a(9)*h2(t) ;
    v(t) = - 0.5 * ( log(2*pi) + log(h1(t)+h2(t)) + (data(t)-mu(t))^2 / (h1(t)+h2(t)) );
end
val = - sum(v); % log likelihood
v = - v; %conditional log likelihoods
end 
