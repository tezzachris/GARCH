function a = thetatoparam(theta)
    

    omega = theta(1); 
    alpha = theta(2); 
    gamma = theta(3); 
    beta = theta(4);
    lambda = theta(5);

    a(1) = exp(omega);
    a(2) = exp(alpha)/(1+exp(alpha));
    a(3) = (gamma);
    a(4) = exp(beta)/(1+exp(beta));
    a(5) = (lambda);

end

