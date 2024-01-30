function theta = paramtotheta(params)
    
    omega = params(1);
    alpha = params(2);
    gamma = params(3);
    beta = params(4);
    lambda = params(5);

    tomega= log(omega);
    talpha= log(alpha/(1-alpha));
    tgamma= (gamma);
    tbeta = log(beta/(1-beta));
    tlambda= (lambda);
    theta = [tomega, talpha, tgamma, tbeta, tlambda]; 

end