function fval = intcalloption(fi,a,b,om,ht1,lam,gam,T,r,K,S)
fip = [i*fi+1; i*fi];

A=fip*r;
B=lam*fip + .5*fip.^2;

% recurse backwards until time t=0
for t=T:-1:1
    Ap=A;
    Bp=B;
    A = Ap + fip*r + Bp*om - .5*log(1-2*a*Bp);
    B = fip*(lam+gam)-.5*gam^2 + b*Bp + (.5*(fip-gam).^2)./(1-2*a*Bp); 
end
    f1=S.^fip(1,:).*exp(A(1,:)+B(1,:)*ht1); 
    f2=S.^fip(2,:).*exp(A(2,:)+B(2,:)*ht1);
    % First integrand
    f1val=real((K.^(-i*fi)).*f1./(i*fi));
    % second integrand
    f2val=real((K.^(-i*fi)).*f2./(i*fi));
    % combined integrand
    fval = (f1val/K - f2val);
end