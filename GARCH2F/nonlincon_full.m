function [c,ceq] = nonlincon_full(a)

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

    B = [b11 + a11*g1^2, b12 + a12*g1^2; 
         b21 + a21*g2^2, b22  + a22*g2^2];

    [~,d]=eig(B); %v=eigenvector, d=eigenvalue

    c(1) = max(abs(diag(d))) - 1;

    ceq = [];
end


