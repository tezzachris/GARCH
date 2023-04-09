
function [VCV] = hessian_2sided(x,data,indexes)
% Computes 2-sided finite difference Hessian

if size(x,2)~=1
    x = x'; %must be column vector
end

n = size(x,1); % x vector 0.03 0.02 0.8 contrained params
estim_flag = 0;
[fx,like] = log_likelihood_chris(x,data,estim_flag,indexes);
t=length(like); %like are the log likelihoods at each point not the sum

% Compute the stepsize (h)
h = eps.^(1/3)*max(abs(x),1e-8);
%xh = x+h;
%h = xh-x; useless step
ee = sparse(diag(h));

% Compute forward and backward steps
gp = zeros(n,1);
gm = zeros(n,1);%like val for me
likep=zeros(t,n); %like v for me
likem=zeros(t,n);
for i=1:n
    [gp(i),likep(:,i)] = log_likelihood_chris(x+ee(:,i),data,estim_flag,indexes);
    [gm(i),likem(:,i)] = log_likelihood_chris(x-ee(:,i),data,estim_flag,indexes);
end

scores=zeros(t,n); 
for i=1:n
    scores(:,i)=(likep(:,i)-likem(:,i))./(2*h(i));
end

hh=h*h';
Hm=NaN*ones(n);
Hp=NaN*ones(n);
% Compute "double" forward and backward steps
for i=1:n
    for j=i:n
        Hp(i,j) = log_likelihood_chris(x+ee(:,i)+ee(:,j),data,estim_flag,indexes);
        Hp(j,i) = Hp(i,j); %symmetric
        Hm(i,j) = log_likelihood_chris(x-ee(:,i)-ee(:,j),data,estim_flag,indexes);
        Hm(j,i) = Hm(i,j);
    end
end
%Compute the hessian
H = zeros(n);
for i=1:n
    for j=i:n
        H(i,j) = (Hp(i,j)-gp(i)-gp(j)+fx+fx-gm(i)-gm(j)+Hm(i,j))/hh(i,j)/2;
        H(j,i) = H(i,j);
    end
end

hess=H/t;
hessinv=hess^(-1);
B=cov(scores);
VCV=(hessinv*B*hessinv)/t;

end 

