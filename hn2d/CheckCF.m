
T=20;
St0=20;
a=[3.32e-5	9.2e-07	25.1	0.95	0.351659507062997];
rf=1e-5;

L=1e6;
logST=zeros(L,1); ht1=logST;
for j = 1:L
    [logSt,ht]=simulate(T+1,a,rf,St0); 
    %Let T=4 so [S0,S1,S2,S3,S4] is T+1 simulated points 
    logST(j)=logSt(end);
    ht1(j)=ht(end);
end

STsim=mean(exp(logST));

phi=1;
h_t0=(a(1) + a(2))/(1-a(2)*a(3)^2-a(4));
h_t1= a(1) + a(4) * h_t0; %is wrong here 

%Let T=4 so [0,1,2,3,4] 
%We know A4,B4 and A3,B3
%We need to go 3->2, 2->1, 1->0 so is T-1 operations

cf=char_fun(phi,T,h_t1,rf,a); 
closedformula=St0*cf(T-1,1);
abs(STsim-closedformula)
