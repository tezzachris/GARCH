clear all

%load results_HN_13_01_2021
load results_HN_31_05_2021
load OptionVet_CALLandPUT_OOM_2020
%load OptionVet_CALLandPUT_OOM_2021

%%% Risk-Neutralized Parameters of Heston-Nandi
%h_HN_1 = 1.761982356863061e-05;
h_HN_1 = theta_joint_mle_HN_GARCH(6);

r_HN = 1.794956716592199e-04;


lambda_HN = -0.5;

omegaHN = theta_joint_mle_HN_GARCH(1);
alphaHN = theta_joint_mle_HN_GARCH(2);
betaHN = theta_joint_mle_HN_GARCH(3);
gammaHN  = theta_joint_mle_HN_GARCH(4) + 0.5 +  theta_joint_mle_HN_GARCH(5);

param_HN_h = [omegaHN, alphaHN, betaHN, gammaHN];

N = 4096*2;

% Real space
x_min = -7.5; x_max = 7.5;
dx=(x_max-x_min)/(N-1);
x=x_min:dx:x_max;

% Fourier space
w_max=pi/dx;
dw=2*w_max/N;
w=[0:dw:w_max, -w_max+dw:dw:-dw];

C_T_HNmod_vet = [];

CFvet = GARCH_FFT_fast_HN_CFvet(1i*w, 999, h_HN_1, r_HN, param_HN_h);

%%% FFT
for j = 1:max(size(C_T_vet))
    S_0 = S_0_vet(j);
    K = K_vet(j);

    Tt = Tt_vet(j);
    dum = OptionType_vet(j);

        % Option payoff
        s = S_0*exp(x);
        
        if dum > 0
            v_call = max(s-K,0);
            % FST method
            CFvet_norm = CFvet(Tt,:);
            v_call = real(ifft(fft(v_call).*CFvet_norm));
            % Interpolate option prices
            C_T_HNmod_vet(j) = interp1(s,v_call,S_0,'PCHIP');
            if C_T_HNmod_vet(j) < 0
                C_T_HNmod_vet(j) = 0;
            end
            % Implied Volatility
            IV_MarketCall(j) = blsimpv(S_0_vet(j), K_vet(j), r_HN*252, Tt/252, C_T_vet(j)...
                                         , 'Limit',1,'Yield',0,'Class', {'Call'});
            IV_HNModCall(j) = blsimpv(S_0_vet(j), K_vet(j), r_HN*252, Tt/252, C_T_HNmod_vet(j)...
                                         , 'Limit',1,'Yield',0,'Class', {'Call'}); 
                                     
%                  cmoney = S_0_vet(j)/K_vet(j);
%                  if cmoney < 0.8 || cmoney > 0.9
%                      IV_MarketCall(j) = 0;
%                      IV_HNModCall(j) = 0;
%                  end
             
                if Tt < 150 || Tt > 365
                    IV_MarketCall(j) = 0;
                    IV_HNModCall(j) = 0;
                end
        else
            v_put = max(K-s,0);
            % FST method
            CFvet_norm = CFvet(Tt,:);
            v_put = real(ifft(fft(v_put).*CFvet_norm));
            % Interpolate option prices
            C_T_HNmod_vet(j) = interp1(s,v_put,S_0,'PCHIP');
            if C_T_HNmod_vet(j) < 0
                C_T_HNmod_vet(j) = 0;
            end
            % Implied Volatility
            IV_MarketPut(j) = blsimpv(S_0_vet(j), K_vet(j), r_HN*252, Tt/252, C_T_vet(j)...
                                         , 'Limit',1,'Yield',0,'Class', {'Put'});
            IV_HNModPut(j) = blsimpv(S_0_vet(j), K_vet(j), r_HN*252, Tt/252, C_T_HNmod_vet(j)...
                                         , 'Limit',1,'Yield',0,'Class', {'Put'});
                                     
%                  pmoney = S_0_vet(j)/K_vet(j);
%                  if pmoney < 0.8 || pmoney > 0.9
%                      IV_MarketPut(j) = 0;
%                      IV_HNModPut(j) = 0;
%                  end
             
                if Tt < 150 || Tt > 365
                     IV_MarketPut(j) = 0;
                     IV_HNModPut(j) = 0;
                end
        end
        
        fprintf("Loop status, j = %d\r\n ", j);
end
    
%%% IMPLIED VOLATILITY RELATIVE ERRORS
RIVCall_HNMod_SE = (nonzeros(IV_MarketCall - IV_HNModCall))./nonzeros(IV_MarketCall);
RIVPut_HNMod_SE = (nonzeros(IV_MarketPut - IV_HNModPut))./nonzeros(IV_MarketPut);
RIVSE_HNMod = [RIVCall_HNMod_SE(:)' RIVPut_HNMod_SE(:)'];

%%% IMPLIED VOLATILITY ROOT MEAN SQUARE ERRORS
IVCall_HNMod_SE =(nonzeros(IV_MarketCall - IV_HNModCall).^2);
IVPut_HNMod_SE = (nonzeros(IV_MarketPut - IV_HNModPut).^2);
IVSE_HNMod = [IVCall_HNMod_SE(:)' IVPut_HNMod_SE(:)'];

IVRMSE_HNMod = sqrt(nanmean(IVSE_HNMod))
   
%%%
% sqrt(mean( ((C_T_vet-C_T_HNmod_vet)./C_T_vet).^2 ))
% sqrt(mean( ((C_T_vet-C_T_HNmod_vet)).^2 ))
% 
% sqrt(mean( (nonzeros(IV_MarketCall-IV_HNModCall)).^2 ))
% sqrt(mean( (nonzeros(IV_MarketPut-IV_HNModPut)).^2 ))
