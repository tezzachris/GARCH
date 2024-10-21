

function [val,v, h1, h2, z1, z2] = loglik_full(theta,data,estim_flag,rf) 

    if estim_flag == 1
        a = thetatoparam_full(theta);
    else
        a = theta;
    end
    
    h1 = zeros(length(data), 1); v=h1; mu = h1; z1=h1; z2=h1; h2 = h1;
    
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

    omega = [w1 ; w2];

    beta = [b11  + a11*g1^2 , b12 + a12*g2^2;
            b21 + a21*g1^2, b22 + a22*g2^2];
     
    I = eye(2,2);

    for t = 1:length(data)
        if t==1
            %Mean as starting value
            h0 = (I-beta) \ (omega + [a11+a12; ...
                                      a22+a21]);
            h1(t) = h0(1);
            h2(t) = h0(2);

        elseif t > 1
                     
            e1 =  ( (data(t-1) - mu(t-1) ) * h1(t-1)) ./( h1(t-1) + h2(t-1) ) ; %et1
            e1sq = e1.^2   ; %et1 squared
            z1(t-1) = e1 ./ sqrt(h1(t-1)); %zt1
            z1sq = e1sq ./ h1(t-1); %zt1 squared
            
            %EKF
            e2 =  data(t-1) - mu(t-1) - e1  ; %et2 
            e2sq = e2.^2   ; %et2 squared
            z2(t-1) = e2 ./ sqrt(h2(t-1)); %zt2   
            z2sq = e2sq ./h2(t-1);

            h1(t) = omega(1) + a11 * ( z1sq - 2*g1*sqrt(h1(t-1))*z1(t-1) + g1^2 * h1(t-1)) ...
                      + b11 * h1(t-1) ...
                      + b12 * h2(t-1) ...
                      + a12 * ( z2sq - 2*g2*sqrt(h2(t-1))*z2(t-1) + g2^2 * h2(t-1) );
        
         
            h2(t) = omega(2) + a22 * ( z2sq - 2*g2*sqrt(h2(t-1))*z2(t-1)  + g2^2 * h2(t-1)) ...
                      + b22 * h2(t-1) ...
                      + b21 * h1(t-1) ...
                      + a21 * ( z1sq - 2*g1*sqrt(h1(t-1))*z1(t-1) + g1^2 * h1(t-1) );
        end
        mu(t) = rf + lam1 * h1(t) + lam2 * h2(t) ;
        v(t) = - 0.5 * ( log(2*pi) + log(h1(t)+h2(t)) + (data(t)-mu(t))^2 / (h1(t)+h2(t)) );
    end
    
    val = - sum(v); % log likelihood
end 
