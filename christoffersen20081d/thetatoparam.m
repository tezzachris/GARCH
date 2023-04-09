function a = thetatoparam(theta)
    
    a(1) = exp(theta(1));
    a(2) = exp(theta(2))/(1+exp(theta(2)));
    a(3) = exp(theta(3));
    a(4) = exp(theta(4))/(1+exp(theta(4)));
    a(5) = exp(theta(5))/(1+exp(theta(5)));
    a(6) = exp(theta(6));
    a(7) = (exp(theta(7))/(1+exp(theta(7))));
    %a(7) = (1-a(6)^2*a(5))*(exp(theta(7))/(1+exp(theta(7))));
    a(8) = (theta(8));

end