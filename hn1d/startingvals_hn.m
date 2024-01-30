function newparams = startingvals_hn(ret,rf)
   
    L = 1000; %multipli di 10
 
    values = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10];

    omega =  rand(L,1) .* values(randi(length(values),[1,L]))'; 
    alpha1 =  rand(L,1) .* values(randi(length(values),[1,L]))';
    beta1 =  0.996*rand(L,1);
    gamma1 =  200 * randn(L,1);
    lambda1 = randn(L,1) ;
 
    a = [omega,alpha1,gamma1,beta1,lambda1]';

    lls = log_likelihood_hn_multi(a,ret,rf,L);
    
    massimo = min( lls((imag(lls)==0)) );
    indice = lls == massimo;
    newparams = a(:,indice);

end


