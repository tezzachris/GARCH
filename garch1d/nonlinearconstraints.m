function [c,ceq] = nonlinearconstraints(a)
            c(1) = a(2) + a(3) - 1  ; %beta1+alpha1*gamma1^2 <= 1
            ceq = [];
end