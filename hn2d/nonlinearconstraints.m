function [c,ceq] = nonlinearconstraints(a)
    c(1) = a(4) + a(2)*a(3)^2  - 1; %beta1+alpha1*gamma1^2 < 1
    c(2) = a(7) + a(5)*a(6)^2 -1; %beta2+alpha2*gamma2^2 < 1
    ceq = [];
end