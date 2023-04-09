function a = thetatoparam(theta)
    
    a(1) = exp(theta(1));
    a(2) = exp(theta(2))/(1+exp(theta(2)));
    a(3) = exp(theta(3));
    a(4) = exp(theta(4))/(1+exp(theta(4)));
    a(5) = (theta(5));

end