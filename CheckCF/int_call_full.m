function fval = int_call_full(fi,T,S,K,h_t1,rf,param)

h1_t = h_t1(1);
h2_t = h_t1(2);

w1 = param(1); a11 = param(2); g1 = param(3); b11 = param(4); 
a22 = param(5); g2 = param(6); b22 = param(7); 
a12 = param(8); a21 = param(9);
b12 = param(10); b21 = param(11);
lam1 = param(12); lam2 = param(13);
w2 = param(14);

phi = [1i*fi+1; 
        1i*fi];

% A and B1, B2 at time T-1
A_old =  phi*rf;
B1_old = phi*lam1 + 0.5*phi.^2;
B2_old = phi*lam2 + 0.5*phi.^2;

% Recurse backwards until time t=0
for t=T-1:-1:1
    A_new = A_old + phi*rf + B1_old*w1 + B2_old*w2  ...
            - 0.5*log(1-2*(a11*B1_old+a21*B2_old)) ...
            - 0.5*log(1-2*(a12*B1_old+a22*B2_old));
    
    B1_new = phi*lam1 ...
            + (b11+a11*g1^2)*B1_old ...
            + (b21+a21*g1^2)*B2_old ...
            + (2*g1^2).*( phi./(2*g1) - (a11*B1_old+a21*B2_old) ).^2 ./(1-2*(a11*B1_old+a21*B2_old));
   
    B2_new = phi*lam2 ...
            + (b12+a12*g2^2)*B1_old ...
            + (b22+a22*g2^2)*B2_old ...
            + (2*g2^2).*( phi./(2*g2) - (a12*B1_old+a22*B2_old) ).^2 ./(1-2*(a12*B1_old+a22*B2_old));
    
    A_old = A_new;
    B1_old = B1_new;
    B2_old = B2_new;

end
    % First CF
    f1=S.^phi(1,:).*exp(A_new(1,:)+B1_new(1,:)*h1_t+B2_new(1,:)*h2_t);
    % Second CF
    f2=S.^phi(2,:).*exp(A_new(2,:)+B1_new(2,:)*h1_t+B2_new(2,:)*h2_t); 

    % First integrand
    f1val=real((K.^(-1i*fi)).*f1./(1i*fi));
    % second integrand
    f2val=real((K.^(-1i*fi)).*f2./(1i*fi));
    % combined integrand
    fval = (f1val - f2val*K); %risk neutral
    
end