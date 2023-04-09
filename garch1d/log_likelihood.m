
function [val,v] = log_likelihood(theta,data,estim_flag,rf) %theta are the transformed params

    if estim_flag == 1
        a = thetatoparam(theta);
    elseif estim_flag == 0
        a = theta;
    end
   
    h = zeros(length(data), 1);
    h(1) = a(1)/(1-a(2)-a(3)) ;
    v=zeros(length(data), 1); mu = v;
    mu(1) = rf + a(4) * h(1);
    v(1)= - 0.5 * ( log(2*pi) + log(h(1)) + (data(1)-mu(1))^2 / h(1) );

    for t = 2:length(data)
        h(t) = a(1) +  a(2) * (data(t-1)-mu(t-1))^2 + a(3) * h(t-1);
        mu(t) = rf + a(4) * h(t);
        v(t) = - 0.5 * (log(2*pi) + log(h(t)) + (data(t)-mu(t))^2 / h(t));
    end
    val = - sum(v); %we want to max so we put minus
    v = - v ; 
end 
