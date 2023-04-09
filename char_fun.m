function [CFvet] = char_fun(phi,T,h_t1,rf,param)

nw = max(size(phi));

CFvet = zeros(T-1, nw); 

omega  = param(1);
alpha = param(2);
gamma = param(3);
beta = param(4);
lambda = param(5);

% A and B and time T-1
A_old = zeros(1,nw);
B_old = zeros(1,nw);

A_old(1,:) = phi*rf;
B_old(1,:) = phi*(lambda+gamma) - 0.5*gamma^2 + 0.5*(phi-gamma).^2;

% Recursion back to time t
for t = 2:T

    A_new = A_old + phi*rf + B_old*omega - 0.5*log(1-2*alpha*B_old);
    B_new = phi*(gamma+lambda) - 0.5*gamma^2 + beta*B_old + (0.5*(phi-gamma).^2)./(1-2*alpha*B_old);
    
    if B_new*alpha > 1/2
        'expectation trick not valid'; 
    end

    A_old = A_new;
    B_old = B_new;
    
    CFvet(t-1,:) = exp(A_new+B_new*h_t1);
end

end