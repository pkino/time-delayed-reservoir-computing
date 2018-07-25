function [NRMSE, NRMSE_C, seedREC] = trialParaTDRC(trial, initNo, paraNum, paraParam, order, gammaCheck, searchCheck)
% paraParam=[theta, learnDimension, biasCheck, inputCheck, a, b, c, p, gamma]
% paraParam(1,:) = [theta, learnDimension, biasCheck, inputCheck, NRMSE(b,1:5)]
if gammaCheck == 1
    
    %  工事中 %
    
    if Model == 'L'
        c = 0;
        gammaInit = 3;
        gammaNum = 3;
        gammaFin = -3;
        gammaData = logspace(gammaInit, gammaFin, gammaNum)';
    elseif Model == 'NL'
        %     c= -0.1;
        gammaInit = -4;
        gammaFin = 1.5;
        %     gammaInter = 0.499;
        %     gammaNum = int32((gammaFin-gammaInit)/gammaInter+1);
        gammaNum = 20;
        %     gammaData = (gammaInit:gammaInter:gammaFin)';
        gammaData = logspace(gammaInit, gammaFin, gammaNum);
    else
        error('モデルを正しく選択してください');
    end
else
    gammaNum = 1;
    if size(paraParam,1) ~= paraNum
        error('パラメータを正しく入力してください');
    end
end

saveParaNum = 3;
seedREC = zeros(paraNum, trial);
NRMSE = NaN(gammaNum, saveParaNum+trial);
NRMSE_C = NaN(gammaNum, saveParaNum+trial);

RCLen = 1500;
dataLen = RCLen + 500;

parfor stepTrial = 1:trial
    %% 入力・目標データの生成
    seed_dataGen = stepTrial+initNo-1;
    [ul, ug, ut, Yl, Yg, Yt] = dataGenerator_NARMA(dataLen,seed_dataGen,order);
    seed_no = stepTrial+initNo; % seed_noを使ったら1足す
    
    %% リザーバ計算
    para_x_kl = cell(paraNum); para_x_kt = cell(paraNum);
    for stepPara = 1:paraNum
        [para_x_kl{stepPara}, para_x_kt{stepPara}] = ...
            timeDelayReservoir(seed_no, ul, ut, paraParam(stepPara,1),paraParam(stepPara,2), ...
            paraParam(stepPara,3), paraParam(stepPara,4), paraParam(stepPara,5), paraParam(stepPara,6), ...
            paraParam(stepPara,7), paraParam(stepPara,8), paraParam(stepPara,9));
        seedREC(stepPara,stepTrial) = seed_no;
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

if searchCheck == 0
    clear ul ug ut Yl Yg Yt;
    
    Date = datestr(datetime('now'),'yyyymmddHHMM');
    save(strcat(Date,'paraTDRC',num2str(paraNum), '_NARMA',num2str(order), '_gammaCheck=', num2str(gammaCheck), '_trial=', num2str(trial), '.mat'), '-v7.3');
end
end

