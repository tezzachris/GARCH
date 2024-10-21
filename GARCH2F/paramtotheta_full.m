function theta = paramtotheta_full(params)
    
    omega1 = params(1);
    omega2 = params(2);
    alpha1 = params(3);
    gamma1 = params(4);
    beta1 = params(5);
    alpha2 = params(6);
    gamma2 = params(7);
    beta2 = params(8);

    alpha12 = params(9);
    alpha21 = params(10);

    beta12 = params(11);
    beta21 = params(12);

    lambda1 = params(13);
    

    tomega1= log(omega1);
    tomega2= log(omega2);
    talpha1= log(alpha1/(1-alpha1));
    tgamma1= (gamma1);
    tbeta1 = log(beta1/(1-beta1));
    talpha2= log(alpha2/(1-alpha2));
    tgamma2= (gamma2);
    tbeta2 = log(beta2/(1-beta2));

    talpha21= log(alpha21/(1-alpha21));
    talpha12= log(alpha12/(1-alpha12));
    
    tbeta12 = log(beta12/(1-beta12));
    tbeta21 = log(beta21/(1-beta21));

    tlambda1 = lambda1;
    
    theta = [tomega1,tomega2, talpha1, tgamma1, tbeta1, talpha2, tgamma2, tbeta2, talpha12, talpha21, tbeta12, tbeta21, tlambda1]; 

end