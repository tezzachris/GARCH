function theta = paramtotheta(params)
    
    omega = params(1);
    alpha = params(2);
    beta = params(3);
    lambda = params(4);

    tomega= log(omega);
    talpha= log(alpha/(1-alpha));
    tbeta = log(beta/(1-beta));
    tlambda= log(lambda);
    theta = [tomega, talpha, tbeta, tlambda]; 

end