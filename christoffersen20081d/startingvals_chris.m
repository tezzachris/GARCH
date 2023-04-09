function newparams = startingvals_chris(ret,rf)
    %keep same size for gamma and alpha and beta
    beta = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    alpha = [1e-13, 1e-12,1e-11, 1e-10,1e-9,1e-8, 1e-7,1e-6,1e-5];
    gamma = [10 : (200-10)/8 :200];
    rip = length(beta);
    a=[repelem(alpha,rip);repmat(gamma,1,rip)];
    a=[repelem(a,1,rip);repelem(beta,size(a,2))]';
    omega = var(ret) * (0.96 - a(:,3) - a(:,1).*a(:,2).^2);
    
    beta2 = beta * rand(1);
    alpha2 = alpha ;
    gamma2 = gamma*rand(1) ;
    rip = length(beta2);
    b=[repelem(alpha2,rip);repmat(gamma2,1,rip)];
    b=[repelem(b,1,rip);repelem(beta2,size(b,2))]';

    maxlambda = 10;
    a = [omega, a, b, maxlambda*rand(size(a,1),1)];

    fval = zeros(size(a,1),1); 
    estim_flag=1;
    for i = 1:size(a,1)
        fval(i)=log_likelihood_chris(paramtotheta(a(i,:)),ret,estim_flag,rf);
    end
    
    massimo = min( fval((imag(fval)==0)) );
    indice=find(fval == massimo);
    newparams= a(indice,:);

end

