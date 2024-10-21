function [ R , h1, h2] = simulate_full(T,a,rf,nsim,R0)
    
    h1 = zeros(T,nsim);  
    h2 = zeros(T,nsim); 
    R =  zeros(T,nsim);
    
    %importa parametri
    w = a(1); a11 = a(2); g1 = a(3); b11 = a(4); 
    a22 = a(5); g2 = a(6); b22 = a(7);
    a12 = a(8); a21 = a(9);
    b12 = a(10); b21 = a(11);
    l1 = a(12); l2 = l1;
    
    omega = [w;
             w];
    beta = [b11 + a11*g1^2 ,  b12 +  a12*g2^2;
            b21 + a21*g1^2 ,  b22 + a22*g2^2];
     
    I = eye(2,2);

    [~,d]=eig(beta); %v=eigenvector, d=eigenvalue
    
    if (max(abs(d),[],'all') > 1)
        error('instabile')
    end
    

    %Come Punto iniziale volatilità usiamo medie
    htminusone = (I-beta) \ ( omega + [a11+a12;a22+a21] ) ;
    h1tminusone = repelem(htminusone(1),1,nsim);
    h2tminusone = repelem(htminusone(2),1,nsim);
    mu = rf + l1 * h1tminusone + l2 * h2tminusone;
    %Usiamo il filtro 
    Z1tminusone = (R0 - mu).*sqrt(h1tminusone) ./ (h1tminusone + h2tminusone);
    Z2tminusone = (R0 - mu).*sqrt(h2tminusone) ./ (h1tminusone + h2tminusone);

    for t = 1:T
            %Equazioni modello
            h1(t,:) = w + a11 * ( Z1tminusone - g1*sqrt(h1tminusone) ).^2 ...
                        + a12 * ( Z2tminusone - g2*sqrt(h2tminusone) ).^2 ...
                        + b11 * h1tminusone ...
                        + b12 * h2tminusone ;

            h2(t,:) = w + a22 * ( Z2tminusone - g2*sqrt(h2tminusone) ).^2  ...
                        + a21 * ( Z1tminusone - g1*sqrt(h1tminusone) ).^2  ...
                        + b21 * h1tminusone ...
                        + b22 * h2tminusone ;
            %Aggiorniamo 
            h1tminusone=h1(t,:);
            h2tminusone=h2(t,:);
            %Simuliamo
            Z1tminusone=randn(1,nsim);
            Z2tminusone=randn(1,nsim);

            mu = rf + l1 * h1(t,:) + l2 * h2(t,:);
            R(t,:) = mu + Z1tminusone.*sqrt(h1(t,:)) + Z2tminusone.*sqrt(h2(t,:)) ;
    end    
    
end 