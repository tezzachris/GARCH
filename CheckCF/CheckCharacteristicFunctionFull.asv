
%% Check by Simulation Characteristic Equation/MGF

cd '/Users/chris/Desktop/unibo/enzo/code/heston2/double/full/options'
%Import parameters
rf=1e-6; % risk free rate per day continuously compounded 
T=60; %days to maturity  
w1 = 1e-10; a11 = 1e-6; g1 = 100; b11 = 0.8; 
a22 = 2e-6; g2 = 120; b22 = 0.6;
a12 = 1e-8; a21 = 1.2e-8;
b12 = 0.003; b21 = 0.005 ;
l1 = 0.01 ; l2 = 0.02;
w2 = 1e-10;
param = [w1,a11,g1,b11,a22,g2, b22, a12, a21, b12, b21, l1, l2, w2]; %paramter vector
omega = [w1 ; w2];

alpha = [a11+a12; a22+a21];
K=84;
g1star = g1 + 0.5 + l1; 
g2star = g2 + 0.5 + l2; 
l1star = -0.5; 
l2star = -0.5;
param_star = [w1,a11,g1star,b11, a22,g2star,b22, a12,a21, b12, b21, l1star, l2star, w2]; %paramter vector
B = [b11  + a11*g1star^2 , b12 + a12*g2star^2;
     b21  + a21*g1star^2 , b22 + a22*g2star^2];
I = eye(2,2);
h_t0= (I-B) \ (omega + alpha);
h_t1=omega+alpha+B*h_t0;
S = 85; % stock price at time 0


%% Characteristic function in Equation 7

%Analytical
phi=3; 
CFvet=CF_HN_F(phi,T,h_t1,rf,param_star); %caratteristica a meno del prezzo S0
fphi_cf = S.^phi.*CFvet(T,:); %caratteristica

%MC Solution 
L = 1e6; %number of simulations
logST=simulate_full(T,param_star,rf,L,S); 
fphi_mc= mean( exp(logST).^(phi) );

error=abs(fphi_mc-fphi_cf)/fphi_mc;


%% Calculate Risk Neutral Double Heston European Call Price 
% Analytical Solution
fun=@(x) int_call_full(x,T,S,K,h_t1,rf,param_star); %equation 11
integr=integral(fun,0,+Inf,'ArrayValued',false); 
call_cf = 0.5*(S-K*exp(-rf*T)) + exp(-rf*T)/pi * integr;

% Monte-Carlo Simulation 
L = 1e6; %number of simulations
[logST]=simulate_full(T,param_star,rf,L,S); 
call_mc = exp(-rf*T) * mean( max( exp(logST) - K ,0) );

error=abs(call_cf - call_mc)/call_mc;

%% Carr Madan

phi = 1; %for call
alpha = 0.75; % used for call option (see [3]), must be positive for put option in [3]
k = log(K);

% FFT with Simpson's rule (as suggested in [2])
N = 2^12; 
eta = 0.25; %small for fine grid
j = 1:N;
v =(j-1)*eta;

lambda = 2*pi/(N*eta); %EQ23 spacing size
b = N*lambda/2; %EQ20
ku = -b+lambda.*(0:N-1); %EQ19
u =  v - (phi.*alpha+1)*1i; %EQ6

CFvet = CF_HN_F( 1i * u, 999, h_t1, rf, param_riskn);
charFunc = S.^(1i * u).*CFvet(T,:);
F = charFunc*exp(-rf*T)./(alpha^2 + phi.*alpha - v.^2 + 1i*(phi.*2*alpha +1).*v); %EQ6
  
SimpsonW = 1/3*(3 + (-1).^(1:N) - [1, zeros(1,N-1)]);
FFTFunc = exp(1i*b*v).*F*eta.*SimpsonW;
payoff = real(fft(FFTFunc));
OptionValue = exp(-phi.*ku.*alpha).*payoff./pi; %EQ24

% Interpolate to get option price for a given strike
y = interp1(ku,OptionValue,k);

error = abs(y - call_mc)/call_mc;

%% Compute Call Option Price via FST

% Griglia
M = 10; %usiamo una discretizzazione con N^M punti spaziali
N=2^M;
L = 7.5; %ampiezza intervallo
xmin=-L; xmax=L;
k=(0:N-1);
dx=(xmax-xmin)/N;
x=xmin+k*dx;
wmin=-pi/dx; wmax=pi/dx; 
dw=2*wmax/N;
w=[0:dw:wmax, wmin+dw:dw:-dw];

%Payoff
s=S*exp(x);
v_call=max(0,s-K); 
v_cap_call=fft(v_call);
CFvet = CF_HN_F(1i*w, 999, h_t1, rf, param_riskn);
CFvet_norm = CFvet(T,:);
v_call=real(ifft(v_cap_call.*CFvet_norm.*exp(-rf*T)));
prezzo_call=interp1(s,v_call,S,'Pchip');

error = abs(prezzo_call - call_mc)/call_mc;
