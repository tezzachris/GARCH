function newparams = start_full(ret,rf)
    
    L = 1000; %multipli di 10

    maxgamma = 600; %can be negative also
    
    values = [1e-6,1e-7,1e-8,1e-9,1e-10];

    omega1 = rand(L,1) .* values(randi(length(values),[1,L]))'; 
    omega2 = rand(L,1) .* values(randi(length(values),[1,L]))'; 
    beta1 =  rand(L,1) ;
    alpha1 = rand(L,1) .* values(randi(length(values),[1,L]))'; 
    gamma1 = rand(L,1) .* maxgamma;
    
    beta2  = rand(L,1);
    alpha2 = rand(L,1) .* values(randi(length(values),[1,L]))'; 
    gamma2 = rand(L,1) .* maxgamma;
  
    alpha21 = rand(L,1) .* values(randi(length(values),[1,L]))'; 
    alpha12 = rand(L,1) .* values(randi(length(values),[1,L]))'; 

    beta12  = 0.5*rand(L,1);
    beta21  = 0.5*rand(L,1);
    
    lambda1 =  5*randn(L,1);
    
    a = [omega1,omega2,alpha1,gamma1,beta1,alpha2,gamma2, beta2, alpha12, alpha21, beta12, beta21, lambda1]';

    lls = loglik_full_multi(a,ret,rf,L);
    
    massimo = min( lls((imag(lls)==0)) );
    indice = lls == massimo;
    newparams = a(:,indice);

end




