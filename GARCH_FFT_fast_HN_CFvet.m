function [CFvet] = GARCH_FFT_fast_HN_CFvet(phi,Tt,h_t1,rF,param)

nw = max(size(phi));

CFvet = zeros(Tt-1, nw);

alpha = param(2);   beta = param(3);   gamma = param(4);
omega  = param(1);

% A and B and time T-1
A_old = zeros(1,nw);
B_old = zeros(1,nw);

A_old(1,:) = phi*rF;
B_old(1,:) = -0.5*phi + 0.5*phi.^2;

% Recursion back to time t
for t = 2:Tt
    A_new = A_old + phi*rF + B_old*omega - 0.5*log(1-2*alpha*B_old);
    B_new = phi*(gamma-0.5) - 0.5*gamma^2 + beta*B_old + (0.5*(phi-gamma).^2)./(1-2*alpha*B_old);
    
    A_old = A_new;
    B_old = B_new;
    
    CFvet(t-1,:) = exp(A_new+B_new*h_t1);
end
end