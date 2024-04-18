function y = HestonFFTVanilla(phi,S,K,T,r,rf,kappa,theta,sigma,rho,v0,alpha,method)
%HESTONFFTVANILLA European FX option price in the Heston model 
%(Carr-Madan approach).
%   Y = HESTONFFTVANILLA(PHI,S,K,T,R,RF,KAPPA,THETA,SIGMA,RHO,V0) returns
%   the price of a European call (PHI=1) or put (PHI=-1) option given 
%   spot price S, strike K, time to maturity (in years) T, domestic R and 
%   foreign RF interest rates, rate of mean reversion KAPPA, average level 
%   of volatility THETA, volatility of volatility SIGMA, correlation 
%   between the Wiener increments driving the spot and vol processes RHO, 
%   and initial volatility VO.

%   Y = HESTONFFTVANILLA(...,ALPHA,METHOD) allows to specify the coefficient 
%   ALPHA of the exponential smoother (default: ALPHA=0.75 for calls and 
%   1+ALPHA=1.75 for puts) and the integration method of the complex
%   integral: METHOD = 0 (default) -> adaptive Gauss-Kronrod quadrature,
%   or METHOD = 1 -> FFT + Simpson's rule (as in [2]).
%    
%   Sample use:
%     >> HestonFFTVanilla(1,100,100,1,.2,.1,.1,.9,.8,0,.05) 
%
%   Reference(s):   
%   [1] H.Albrecher, P.Mayer, W.Schoutens, J.Tistaert (2006) The little
%       Heston trap, Wilmott Magazine, January: 83–92.
%   [2] P.Carr, D.Madan (1998) Option valuation using the Fast Fourier 
%       transform, J. Computational Finance 2, 61-73.
%   [3] A.Janek, T.Kluge, R.Weron, U.Wystup (2010) FX smile in the
%       Heston model, see http://ideas.repec.org/p/pra/mprapa/25491.html
%       {Chapter prepared for the 2nd edition of "Statistical Tools for 
%       Finance and Insurance", P.Cizek, W.Härdle, R.Weron (eds.), 
%       Springer-Verlag, forthcoming in 2011.}  
%   [4] M.Schmelzle (2010) Option Pricing Formulae using Fourier Transform: 
%       Theory and Application, Working paper.
   
%   Written by Agnieszka Janek (2010.07.23)
%   Revised by Rafal Weron (2010.10.08)
%   Revised by Agnieszka Janek and Rafal Weron (2010.10.21, 2010.12.27)

if (nargin < 11)
  error ('Wrong number of input arguments.')
else
  if (nargin < 13)
    % Set default value of method: 
    method = 0; % adaptive Gauss-Kronrod quadrature
    if (nargin == 11)
      % Set parameter alpha (see [3])
      alpha = 0.75; % used for call option
    end
  end
end

if (phi==-1) %put option
  alpha = alpha + 1; 
end

s0 = log(S);
k = log(K);
    
if (method == 0)
  % Integrate using adaptive Gauss-Kronrod quadrature
  y = exp(-phi.*k.*alpha).*quadgk(@(v) HestonFFTVanillaInt(phi,s0,k,T,r,rf,kappa,theta,sigma,rho,v0,alpha,v),0,inf,'RelTol',1e-8)./pi;
else


  % FFT with Simpson's rule (as suggested in [2])
  N = 2^10; 
  eta = 0.25;
  v =(0:N-1)*eta;
    
  lambda = 2*pi/(N*eta);
  b = N*lambda/2;
  ku = -b+lambda.*(0:N-1);
  u = v - (phi.*alpha+1)*1i;

  d = sqrt((rho.*sigma.*u*1i-kappa).^2+sigma.^2.*(1i*u+u.^2));
  g = (kappa-rho.*sigma.*1i.*u-d)./(kappa-rho.*sigma*1i.*u+d);
    
  % Characteristic function (see [1])
  A = 1i*u.*(s0 + (r-rf).*T);
  B = theta.*kappa.*sigma.^(-2).*((kappa-rho.*sigma*1i.*u-d).*T-2*log((1-g.*exp(-d.*T))./(1-g)));
  C = v0.*sigma.^(-2).*(kappa-rho.*sigma*1i.*u-d).*(1-exp(-d.*T))./(1-g.*exp(-d.*T));
  charFunc = exp(A + B + C); 
  F = charFunc*exp(-r*T)./(alpha^2 + phi.*alpha - v.^2 + 1i*(phi.*2*alpha +1).*v);
    
  % Use Simpson's approximation to calculate FFT (see [2])
  SimpsonW = 1/3*(3 + (-1).^[1:N] - [1, zeros(1,N-1)]);
  FFTFunc = exp(1i*b*v).*F*eta.*SimpsonW;
  payoff = real(fft(FFTFunc));
  
  %EQUATION 3
  OptionValue = exp(-phi.*ku.*alpha).*payoff./pi;
  % Interpolate to get option price for a given strike
  y = interp1(ku,OptionValue,k);
end

%%%%%%%%%%%%% INTERNALLY USED ROUTINE %%%%%%%%%%%%%

function payoff = HestonFFTVanillaInt(phi,s0,k,T,r,rf,kappa,theta,sigma,rho,v0,alpha,v)
%HESTONFFTVANILLAINT Auxiliary function used by HESTONFFTVANILLA.
%   PAYOFF=HESTONFFTVANILLAINT(phi,s0,k,T,r,rf,kappa,theta,sigma,rho,v0,alpha,v)
%   returns the values of the auxiliary function evaluated at points V, 
%   given log(spot price) S0, log(strike) K, time to maturity (in years) T,
%   domestic interest rate R, foreign interest rate RF, level of mean
%   reversion KAPPA, long-run variance THETA, vol of vol SIGMA, correlation 
%   RHO, initial volatility VO, damping coefficient ALPHA and option type PHI:
%   PHI = 1 --> call option,
%   PHI = -1 --> put option.
   
u = v - (phi.*alpha+1)*1i;
d = sqrt((rho.*sigma.*u*1i-kappa).^2+sigma.^2.*(1i*u+u.^2));
g = (kappa-rho.*sigma.*1i.*u-d)./(kappa-rho.*sigma*1i.*u+d);

% Characteristic function (see [1])
A = 1i*u.*(s0 + (r-rf).*T);
B = theta.*kappa.*sigma.^(-2).*((kappa-rho.*sigma*1i.*u-d).*T-2*log((1-g.*exp(-d.*T))./(1-g)));
C = v0.*sigma.^(-2).*(kappa-rho.*sigma*1i.*u-d).*(1-exp(-d.*T))./(1-g.*exp(-d.*T));
charFunc = exp(A + B + C); 

FFTFunc = charFunc*exp(-r*T)./(alpha^2 + phi.*alpha - v.^2 + 1i*(phi.*2*alpha +1).*v);
payoff = real(exp(-1i.*v.*k).*FFTFunc);

