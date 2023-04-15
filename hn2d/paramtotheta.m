function theta = paramtotheta(params)
    
    omega = params(1);
    alpha1 = params(2);
    gamma1 = params(3);
    beta1 = params(4);
    alpha2 = params(5);
    gamma2 = params(6);
    beta2 = params(7);
    lambda1 = params(8);
    lambda2 = params(9);

    tomega= log(omega);
    talpha1= log(alpha1/(1-alpha1));
    tgamma1= log(gamma1);
    tbeta1 = log(beta1/(1-beta1));
    talpha2= log(alpha2/(1-alpha2));
    tgamma2= log(gamma2);
    tbeta2 = log(beta2/(1-beta2));
    tlambda1= lambda1;
    tlambda2= lambda2; %can be negative
    theta = [tomega, talpha1, tgamma1, tbeta1, talpha2, tgamma2, tbeta2, tlambda1, tlambda2]; 

end