
%Characteristic function analytical solution via Equation 7 of Heston-Nandi

function [CFvet] = CF_HN_Full(phi,T,rf,a)

nw = max(size(phi)); 

CFvet = zeros(T-1, nw); 

    w1 = a(1); 
    w2 = a(2);
    a11 = a(3); 
    g1 = a(4); 
    b11 = a(5); 
    a22 = a(6); 
    g2 = a(7); 
    b22 = a(8); 
    a12 = a(9); 
    a21 = a(10);
    b12 = a(11); 
    b21 = a(12);
    lam1 = a(13);
    lam2 = lam1;


g1 = g1 + 0.5 + lam1; %is g1star
g2 = g2 + 0.5 + lam2; %is g2star

lam1 = -0.5; %is lambdastar1
lam2 = -0.5; %is lambdastar2, change code down when you have time 

omega = [w1 ; w2];
B = [b11 + a11*g1^2 , b12+a12*g2^2;
     b21 + a21*g1^2 , b22 + a22*g2^2];
alpha = [a11+a12; a22+a21];
I = eye(2,2);
h_t0= (I-B) \ (omega + alpha);
h_t1=omega+alpha+B*h_t0;

h1_t = h_t1(1);
h2_t = h_t1(2);

% A and B coefficients at time T
%A_old = zeros(1,nw); B_old = zeros(1,nw);

% A and B coefficients at time T-1
A_old = phi*rf;
B1_old = phi*lam1 + 0.5*phi.^2;
B2_old = phi*lam2 + 0.5*phi.^2;

% Recursion back to time t, 
% from T-1 to 0 we need T-1 steps
for t = 1:T-1

    A_new = A_old + phi.*rf + B1_old*w1 + B2_old*w2  ...
            - 0.5*log(1-2*(a11*B1_old+a21*B2_old)) ...
            - 0.5*log(1-2*(a12*B1_old+a22*B2_old));
    B1_new = phi.*lam1 ...
            + (b11+a11*g1^2)*B1_old ...
            + (b21+a21*g1^2)*B2_old ...
            + (2*g1^2)*( phi./(2*g1) - (a11*B1_old+a21*B2_old) ).^2 ./(1-2*(a11*B1_old+a21*B2_old));
    B2_new = phi.*lam2 ...
            + (b12+a12*g2^2)*B1_old ...
            + (b22+a22*g2^2)*B2_old ...
            + (2*g2^2)*( phi./(2*g2) - (a12*B1_old+a22*B2_old) ).^2 ./(1-2*(a12*B1_old+a22*B2_old));
    
    A_old = A_new;
    B1_old = B1_new;
    B2_old = B2_new;
    
    CFvet(t,:) = exp(A_new+B1_new.*h1_t+B2_new.*h2_t);
end

end