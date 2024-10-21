
function [lls, h1, h2, z1, z2] = loglik_full_multi(a,data,rf,L) 
    
    h1 = zeros(length(data), L); v=h1; mu = h1; z1=h1; z2=h1; h2 = h1;
    
    %Parameters
    w1 = a(1,:); 
    w2 = a(2,:);
    a11 = a(3,:); 
    g1 = a(4,:); 
    b11 = a(5,:); 
    a22 = a(6,:); 
    g2 = a(7,:); 
    b22 = a(8,:); 
    a12 = a(9,:); 
    a21 = a(10,:);
    b12 = a(11,:); 
    b21 = a(12,:);
    l1 = a(13,:); 
    l2 = l1;

    omega = [w1 , w2];

    for t = 1:length(data)
        if t==1
                p1 = (omega(1) + a11 + a12)./(1-b11-a11.*g1.^2) ; 
                p2 = (b12 + a12 .* g2.^2)./(1-b11-a11.*g1.^2) ;
            
                num = omega(2) + (b21 + a21 .* g1.^2) .* p1 + a21 + a22 ;
                den = 1 - p2 .* (a21 .* g1.^2 + b21) - a22.*g2.^2 - b22 ;
                Eh2 = num./den;
            
                Eh1 = p1 + p2 .* Eh2;

                h1(t,:) = Eh1; 
                h2(t,:)=  Eh2; 
                
        elseif t > 1
                
                e1 =  ( (data(t-1) - mu(t-1,:) ) .* h1(t-1,:)) ./( h1(t-1,:) + h2(t-1,:) ) ; 
                e1sq = e1.^2 ;
                z1(t-1,:) =  e1 ./sqrt(h1(t-1,:));
                z1sq = e1sq ./ h1(t-1,:);

                e2 =  data(t-1) - mu(t-1,:) - e1;  
                e2sq = e2.^2  ;
                z2(t-1,:) = e2 ./sqrt(h2(t-1,:));
                z2sq = e2sq ./ h2(t-1,:);
                            
                h1(t,:) = w1 + a11 .* ( z1sq - 2*g1.*sqrt(h1(t-1,:)).*z1(t-1,:) + g1.^2 .* h1(t-1,:)) ...
                            + b11 .* h1(t-1,:) + b12 .* h2(t-1,:) ...
                            + a12 .* ( z2sq - 2*g2.*sqrt(h2(t-1,:)).*z2(t-1,:) + g2.^2 .* h2(t-1,:) );
 
                h2(t,:) = w2 + a22 .* ( z2sq - 2*g2.*sqrt(h2(t-1,:)).*z2(t-1,:) + g2.^2 .* h2(t-1,:))  ...
                            + b22 .* h2(t-1,:) + b21 .* h1(t-1,:) ...
                            + a21 .* ( z1sq - 2*g1.*sqrt(h1(t-1,:)).*z1(t-1,:) + g1.^2 .* h1(t-1,:) );
        end
        mu(t,:) = rf + l1.*h1(t,:) + l2.*h2(t,:) ;
        v(t,:) = - 0.5 * ( log(2*pi) + log(h1(t,:)+h2(t,:)) + (data(t,:)-mu(t,:)).^2 ./ (h1(t,:)+h2(t,:)) );
    end
    lls = - sum(v,1); % log likelihoods
end 
