
function [val, h, v] = log_likelihood_hn_multi(a,data,rf,L) %data is zeta N(0,1) vector Tx1

    h = zeros(length(data), L); v=h; mu = h; z=h;
    
    h(1,:)= (a(1,:)+a(2,:)) ./ (1 - a(2,:).*a(3,:).^2 - a(4,:));
    
    mu(1,:) = rf  + a(5)*h(1,:);

    v(1,:) = - 0.5 * ( log(2*pi) + log(h(1,:)) + (data(1,:)-mu(1,:)).^2 ./ h(1,:) );
    
    for t = 2:length(data)
        z(t-1,:) = (data(t-1,:)-mu(t-1,:) )./sqrt(h(t-1,:));
        h(t,:) =   a(1,:) + a(2,:) .* ( z(t-1,:) - a(3,:).*sqrt(h(t-1,:)) ).^2 + a(4,:) .* h(t-1,:);   
       
        mu(t,:) = rf + a(5,:).*h(t,:) ;
 
        v(t,:) = - 0.5 * ( log(2*pi) + log(h(t,:)) + (data(t,:)-mu(t,:)).^2 ./ h(t,:) );
        
    end
    val = - sum(v,1); % log likelihood
   
end 
