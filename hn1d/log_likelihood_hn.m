
function [val, v, h] = log_likelihood_hn(theta,data,estim_flag,rf) %data is zeta N(0,1) vector Tx1

    if estim_flag == 1
        a = thetatoparam(theta);
    elseif estim_flag == 0
        a = theta;
    end

    h = zeros(length(data), 1); v=h; mu = h; z=h;
    w = a(1); a1 = a(2); g1 = a(3); b1 = a(4); l1 = a(5);
    h(1)= (w+a1)/(1 - a1*g1^2 - b1);
    mu(1) = rf  + l1*h(1);

    v(1) = - 0.5 * ( log(2*pi) + log(h(1)) + (data(1)-mu(1))^2 / h(1) );
    
    for t = 2:length(data)
        z(t-1) = (data(t-1)-mu(t-1) )/sqrt(h(t-1));
        h(t) =   w + a1 * ( z(t-1) - g1*sqrt(h(t-1)) )^2 + b1 * h(t-1);         
        mu(t) = rf + l1*h(t) ;
 
        v(t) = - 0.5 * ( log(2*pi) + log(h(t)) + (data(t)-mu(t))^2 / h(t) );
    end
    val = - sum(v); % log likelihood
   
end 
