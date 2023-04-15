function a = thetatoparam(theta)
    
    a(1) = exp(theta(1)); %omega
    a(2) = exp(theta(2))/(1+exp(theta(2))); %alpha1
    a(3) = exp(theta(3)); %gamma1
    a(4) = exp(theta(4))/(1+exp(theta(4))); %beta1
    a(5) = exp(theta(5))/(1+exp(theta(5))); %alpha2
    a(6) = exp(theta(6)); %gamma2
    a(7) = exp(theta(7))/(1+exp(theta(7))); %beta2
    a(8) = (theta(8)); %lambda1
    a(9) = (theta(9)); %lambda2

end