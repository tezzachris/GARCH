

%Option log-likelihood

%df is dataframe containing option prices, underlying , days to maturity, flag put/call, strike

function [val, v, sigma2u, IVRMSE] = loglik_option_full_iv(theta,df,estim_flag,rf) %data is zeta N(0,1) vector Tx1

    if estim_flag == 1
        a = thetatoparam(theta);
    elseif estim_flag == 0
        a = theta;
    end

N = 4096*2; %grid width

% Real space
x_min = -7.5; x_max = 7.5; % min and max daily returns grid
dx=(x_max-x_min)/(N-1);
x=x_min:dx:x_max;

% Fourier space
w_max=pi/dx;
dw=2*w_max/N;
w=[0:dw:w_max, -w_max+dw:dw:-dw];

CFvet = CF_HN_Full(1i*w, 999, rf, a);
n = size(df,1);
CTmod_vet = zeros(n,1); %option prices for model
IV_HNModCall = zeros(n,1); %implied volatility for model

for j = 1:n

    S = df(j,:).SPX;
    K = df(j,:).Strike;
    Tt = df(j,:).DTM;
    dum = df(j,:).FlagCallPut;
    s = S * exp(x); % Option payoff

    %FST method
    CFvet_norm = CFvet(Tt-1,:);

    if dum == 1
            v_call = max(s-K,0);
            v_call = real(ifft(fft(v_call).*CFvet_norm));
            
            % Interpolate option prices
            CT = interp1(s,v_call,S);
            CTmod_vet(j) = CT * exp(-rf*Tt);

            if CTmod_vet(j) < 0
                CTmod_vet(j) = 0;
            end
            
            IV_HNModCall(j) = blsimpv(S , K , rf * 252, Tt/252, CTmod_vet(j)...
                                        , 'Limit',1,'Yield',0,'Class', {'Call'}); %requires annualized params 

       elseif dum == 0

            v_put = max(K-s,0);
            v_put = real(ifft(fft(v_put).*CFvet_norm));

            %Interpolate option prices
            CT = interp1(s,v_put,S);
            CTmod_vet(j) = CT * exp(-rf*Tt);

            if CTmod_vet(j) < 0
                CTmod_vet(j) = 0;
            end

            IV_HNModCall(j) = blsimpv(S , K , rf * 252, Tt/252, CTmod_vet(j)...
                                        , 'Limit',1,'Yield',0,'Class', {'Put'}); 
                                  
    end

end

idx = find(CTmod_vet~=0);
UiCall = (df.IVMarket(idx) - IV_HNModCall(idx))./df.IVMarket(idx);
sigma2u = var(UiCall,'omitnan');
v = - 0.5 * ( log(2*pi*sigma2u) + (UiCall).^2 ./ sigma2u );

val = - sum(v,'omitnan'); %negative log-likelihood

IVRMSE = 100*sqrt(mean( (df.IVMarket(idx) - IV_HNModCall(idx)).^2, 'omitnan' ));

end 
