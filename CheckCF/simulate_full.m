function [ logSt, h1, h2, z1, z2] = simulate_full(T,a,rf,nsim,St0)
      
    T = T+1; %maturity
    %in this function we save only values at maturity for speed!
    w1 = a(1); 
    a11 = a(2); 
    g1 = a(3); 
    b11 = a(4); 
    a22 = a(5); 
    g2 = a(6); 
    b22 = a(7);
    a12 = a(8); 
    a21 = a(9);
    b12 = a(10); 
    b21 = a(11);
    l1 = a(12); 
    l2 = a(13);
    w2 = a(14);

    omega = [w1;
             w2];

    beta = [b11 + a11*g1^2 ,  b12 +  a12*g2^2;
            b21 + a21*g1^2 ,  b22 +  a22*g2^2];
     
    I = eye(2,2);
    alpha = [a11+a12;a22+a21];
    [~,d]=eig(beta); %v=eigenvector, d=eigenvalue
    if (max(abs(d),[],'all') > 1)
        error('non-stationary')
    end
    
    z1 = randn(T,nsim);  
    z2 = randn(T,nsim);

    %we should start at t=0 
    for t = 1:T
        if t==1
            h0 = (I-beta) \ (omega + alpha);
            h1 = h0(1);
            h2 = h0(2);
            logSt =  log(St0) ;
        else
            h1 = omega(1) + a11 * ( z1(t-1,:) - g1*sqrt(h1) ).^2 ...
                        + a12 * ( z2(t-1,:) - g2*sqrt(h2) ).^2 ...
                        + b11 * h1 ...
                        + b12 * h2 ;
            h2 = omega(2) + a22 * ( z2(t-1,:) - g2*sqrt(h2) ).^2  ...
                        + b22 * h2 ...
                        + b21 * h1 ...
                        + a21 * ( z1(t-1,:) - g1*sqrt(h1) ).^2 ;
            logSt = logSt + rf + l1 * h1  + l2 * h2 ...
                        + sqrt(h1).*z1(t,:) + sqrt(h2).*z2(t,:) ;
        end
        
    end

end 