function [bestNRMSE, NRMSE, NRMSE_C] = ESNTypeTDRC(trial, paraNum, theta, ... 
    oneDim, biasCheck, inputCheck, a, gamma, W,  order,  searchCheck, ang)
% paraParam=[theta, learnDimension, biasCheck, inputCheck, a, b, c, p, gamma]
% paraParam(1,:) = [theta, learnDimension, biasCheck, inputCheck, NRMSE(b,1:5)]

saveParaNum = 3;
NRMSE = NaN(gammaNum, saveParaNum+trial);
NRMSE_C = NaN(gammaNum, saveParaNum+trial);

RCLen = 1500;
dataLen = RCLen + 500;

for stepTrial = 1:trial
    %% 入力・目標データの生成
    [ul, ug, ut, Yl, Yg, Yt] = dataGenerator_NARMA(dataLen,order);
    
    %% リザーバ計算
    para_x_kl = cell(paraNum); para_x_kt = cell(paraNum);
    for stepPara = 1:paraNum
        [para_x_kl{stepPara}, para_x_kt{stepPara}] = ...
            ESNTypeTimeDelayReservoir(stepTrial, ul, ut,theta, oneDim, biasCheck, inputCheck, a, gamma, W);
    end
    
    %% 学習とテスト
    try
        [NRMSE(saveParaNum+stepTrial), NRMSE_C(saveParaNum+stepTrial)] = TDRC(RCLen, vertcat(para_x_kl{:}), vertcat(para_x_kt{:}), Yl, Yt);
    catch
        NRMSE(saveParaNum+stepTrial) =  NaN; NRMSE_C(saveParaNum+stepTrial) =  NaN;
    end
end

NRMSE(1,1) = mean(NRMSE(1,saveParaNum+1:end),'omitnan');
NRMSE(1,2) = std(NRMSE(1,saveParaNum+1:end),0,2);
NRMSE_C(1,1) = mean(NRMSE_C(1,saveParaNum+1:end),'omitnan');
NRMSE_C(1,2) = std(NRMSE_C(1,saveParaNum+1:end),0,2);

[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,1));
bestNRMSE_C_ofBestNRMSE = NRMSE_C(bestNRMSEIndex);
[bestNRMSE_C, bestNRMSE_CIndex] = min(NRMSE_C(:,1));

if searchCheck == 0
    clear ul ug ut Yl Yg Yt;
    
    Date = datestr(datetime('now'),'yyyymmddHHMM');
    
    save(strcat(Date, 'ang', num2str(ang), '_bridge_TDRC',num2str(paraNum), '_NARMA',num2str(order), '_gammaCheck=', num2str(gammaCheck), '_trial=', num2str(trial), '.mat'), '-v7.3');
end
end

