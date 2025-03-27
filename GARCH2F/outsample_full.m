

function [Res,QL] = outsample_full(primov, ultimov, ret, rf, signif )

    nsim = 100000;
   
    subret = ret(primov:ultimov,1);
    n = ultimov - primov + 1;
    ninsample = ceil(0.6 * n);   
    j = 0;
    prevparam = [];
    VARs1 = zeros(n-ninsample,1);  %storing the value at risk
    VARs2=VARs1; VARs3=VARs1; VARs4=VARs1; VARs5=VARs1; 
    
    for i = ninsample : 1 : n-1
        j = j + 1 ; 
        insample = subret(j:i,1);
        [param] = fmincon_full(insample,rf,prevparam);

        prevparam = param;
        R0 = insample(end); 
        simuls = simulate_full(5,param,rf,nsim,R0); %usavo 2 prima meglio 1?

        VARs1(j) = quantile( simuls(1,:), 1-signif);
        VARs2(j) = quantile( simuls(2,:), 1-signif);
        VARs3(j) = quantile( simuls(3,:), 1-signif);
        VARs4(j) = quantile( simuls(4,:), 1-signif);
        VARs5(j) = quantile( simuls(5,:), 1-signif);
        
    end
    
  %Outofsample observations
yt = subret(ninsample+1:n);

   %VaR Back Test
vbt1=varbacktest(  yt , -VARs1, 'VaRLevel',signif); %Note: -VaR and 0.95 signif to test the 0.05 quantile
vbt2=varbacktest(  yt , -VARs2, 'VaRLevel',signif);
vbt3=varbacktest(  yt , -VARs3, 'VaRLevel',signif);
vbt4=varbacktest(  yt , -VARs4, 'VaRLevel',signif);
vbt5=varbacktest(  yt , -VARs5, 'VaRLevel',signif);

%store results
p1 = cc(vbt1);p2 = cc(vbt2);p3 = cc(vbt3);p4 = cc(vbt4);p5 = cc(vbt5);
ResArr = [p1.LRatioPOF,p1.PValuePOF,p1.LRatioCC,p1.PValueCC;
          p2.LRatioPOF,p2.PValuePOF,p2.LRatioCC,p2.PValueCC;
          p3.LRatioPOF,p3.PValuePOF,p3.LRatioCC,p3.PValueCC;
          p4.LRatioPOF,p4.PValuePOF,p4.LRatioCC,p4.PValueCC;
          p5.LRatioPOF,p5.PValuePOF,p5.LRatioCC,p5.PValueCC];

Res = array2table(ResArr,'RowNames',{'1day' '2day' '3day' '4day' '5day'}, ...
                  'VariableNames', {'LRatioPOF' 'PValuePOF' 'LRatioCC' 'PValueCC'});


%Quantile Loss
QL1 =  ( (yt - VARs1) .* (1-signif - (yt < VARs1) ) );
QL2 =  ( (yt - VARs2) .* (1-signif - (yt < VARs2) ) );
QL3 =  ( (yt - VARs3) .* (1-signif - (yt < VARs3) ) );
QL4 =  ( (yt - VARs4) .* (1-signif - (yt < VARs4) ) );
QL5 =  ( (yt - VARs5) .* (1-signif - (yt < VARs5) ) );

QL = [QL1,QL2,QL3,QL4,QL5];


end

