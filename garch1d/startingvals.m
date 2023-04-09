function newparams = startingvals(ret,rf)
   
    beta = [0.3 0.6 0.9];
    alpha = [1e-9 1e-6 1e-3];
    nparams=4;
    maxlambda = 6;
    fval = zeros(length(beta)*length(alpha),1); 
    np=zeros(length(beta)*length(alpha), nparams) ;
    j=0; estim_flag=1;
    for i = beta
            for k = alpha
                j=j+1;
                param=[var(ret)*(0.996 - i - k) , k, i, rand(1)*maxlambda];
                [np(j,:), fval(j)]=fminsearch('log_likelihood' , paramtotheta(param) , optimset('display','off','Maxiter',1000),ret,estim_flag,rf);
            end
    end
    massimo = min( fval((imag(fval)==0)) );
    indice=find(fval == massimo);
    newparams= np(indice,:);

    
%     
%     nmin = 11; nmax = 7;
%     col = zeros(nmin-nmax,10);
%     j=0;
%     for i = nmin : -1 : nmax %1e-nmin -> 1e-nmax
%         j = j + 1;
%         col(j,:) = 1/10^(i+1) : 1/10^(i+1) : 1/10^(i);
%     end
%     col=nonzeros(col); 
%     idx = randi( [1, size(col,1)] , 1, nsim);
%     alpha= col(idx);
% 
%     lambda = rand(nsim,1)*15;
% 
%     for i = 1:nsim
%         ub = .99996; %upperbound
%         lb = 0.1;
%         beta(i) = (ub - lb)*rand(1) + lb;
%         omega(i) = var(data) * (1-beta(i)-alpha(i)) ;
%         params= [omega(i), alpha(i), beta(i), lambda(i)];
%         loglik(i) = log_likelihood(params,data,0,rf);
%     end
%     
%     [~,index]=sort(-loglik,'descend');
%     index=index(1);
%     startparams = [omega(index), alpha(index), beta(index), lambda(index)];

end