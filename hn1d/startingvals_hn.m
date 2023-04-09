function newparams = startingvals_hn(ret,rf)
    
    beta = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    alpha = [1e-16,1e-14,1e-13,1e-12,1e-11,1e-10,1e-8,1e-6,1e-5];
    gamma = [10 30 50 100 150 300 500 750 1000];
    rip = length(beta);
    a=[repelem(alpha,rip);repmat(gamma,1,rip)];
    a=[repelem(a,1,rip);repelem(beta,size(a,2))]';
    omega = var(ret) * (0.6 - a(:,3) - a(:,1).*a(:,2).^2);

    maxlambda = 10; %max absolute max for lambda
    a1 = [omega, a, maxlambda*rand(size(a,1),1) ];
   
    fval = zeros(size(a1,1),1); 
    estim_flag=1;

    for i = 1:size(a1,1)
        %[np(j,:), fval(j)]=fminsearch('log_likelihood' , paramtotheta(param) , optimset('display','off','Maxiter',1000),ret,estim_flag,rf);
        fval(i)=log_likelihood_hn(paramtotheta(a1(i,:)),ret,estim_flag,rf);
    end
    
    massimo = min( fval((imag(fval)==0)) );
    indice=find(fval == massimo);
    newparams= a1(indice,:);

%     %Try also negative lambda's
%     a2 = [omega, a, -1e-5 * rand(size(a,1),1)];
%     fval2 = zeros(size(a2,1),1); 
% 
%     for i = 1:size(a2,1)
%         %[np(j,:), fval(j)]=fminsearch('log_likelihood' , paramtotheta(param) , optimset('display','off','Maxiter',1000),ret,estim_flag,rf);
%         fval2(i)=log_likelihood_hn(paramtotheta(a2(i,:)),ret,estim_flag,rf);
%     end
%     
%     massimo2 = min( fval2((imag(fval2)==0)) );
%     if massimo > massimo2
%         indice=find(fval == massimo);
%         newparams= a1(indice,:);
%     else
%         indice=find(fval2 == massimo2);
%         newparams= a2(indice,:);
%     end
%  


end