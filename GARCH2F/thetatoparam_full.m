function a = thetatoparam_full(theta)
    
    a(1) = exp(theta(1)); %omega
    a(2) = exp(theta(2)); %omega
    a(3) = exp(theta(3))/(1+exp(theta(3))); %alpha1
    a(4) = (theta(4)); %gamma1
    a(5) = exp(theta(5))/(1+exp(theta(5))); %beta1

    a(6) = exp(theta(6))/(1+exp(theta(6))); %alpha2
    a(7) = (theta(7)); %gamma2
    a(8) = exp(theta(8))/(1+exp(theta(8))); %beta2

    a(9) = exp(theta(9))/(1+exp(theta(9))); %alpha12
    a(10) = exp(theta(10))/(1+exp(theta(10))); %alpha21

    a(11) = exp(theta(11))/(1+exp(theta(11))); %beta12
    a(12) = exp(theta(12))/(1+exp(theta(12))); %beta21

    a(13) = theta(13); %lambda1


end