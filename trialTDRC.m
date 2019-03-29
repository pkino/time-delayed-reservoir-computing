function trialTDRC(trial, Model, Task, order, theta, learnDimension, biasCheck, inputCheck, ...
    cMin, cMax, pMin)
dbstop if error

%% モデルの選択
if Model == 'L'
    c = 0;
    gammaInit = 3;
    gammaNum = 3;
    gammaFin = -3;
    gammaData = logspace(gammaInit, gammaFin, gammaNum)';
elseif Model == 'NL'
    %     c= -0.1;
    gammaInit = -2;
    gammaFin = 1;
    %     gammaInter = 0.499;
    %     gammaNum = int32((gammaFin-gammaInit)/gammaInter+1);
    gammaNum = 20;
    %     gammaData = (gammaInit:gammaInter:gammaFin)';
    gammaData = logspace(gammaInit, gammaFin, gammaNum);
else
    error('モデルを正しく選択してください');
end

eigMin=-1.5; eigMax=-6; eigNum=5; gapMax=-2; gapMin=-5.5; gapNum=5; % eigMin=-1.5; eigMax=-9; gapMax=-2; gapMin=-8;
[a_mat, b_mat] = find_ab(eigMin, eigMax, eigNum, gapMax, gapMin, gapNum);

cNum = 20;
% cMin = -0.1;
% cMax = 1.5;
cData = logspace(cMin,cMax,cNum);

p=3;

singleSearchNum = 4;
saveNum = 2;
singleSearchParams = gapNum*eigNum*cNum*pNum*gammaNum;
singleParamsSet = NaN(singleSearchParams,singleSearchNum+3+1);

RCLen = 1500;
dataLen = RCLen + 500;
dataGen = str2func(strcat('dataGenerator_',Task));

index = 1;
for step_gap = 1:gapNum
    for step_eig = 1:eigNum
        for step_c = 1:cNum
            for step_p = 1:pNum
                for step_g = 1:gammaNum
                    singleParamsSet(index, 1:singleSearchNum) = [a_mat(step_gap,step_eig), b_mat(step_gap,step_eig), ...
                        cData(step_c), pData(step_p), gammaData(step_g)];
                    index = index + 1;
                end
            end
        end
    end
end

parfor step = 1:singleSearchParams
    for stepTrial = 1:trial
        %% 入力・目標データの生成
        [ul, ~, ut, Yl, ~, Yt] = dataGen(dataLen, order);
        
        %% リザーバ計算
        [x_kl, x_kt] ...
            = timeDelayReservoir(stepTrial, ul, ut, theta, learnDimension, biasCheck, inputCheck, ...
            singleParamsSet(step,1),singleParamsSet(step,2), singleParamsSet(step,3), p, singleParamsSet(step,4));
        
        %% 学習とテスト
        [N,NC,Ylp,Ytp]  = TDRC(RCLen, x_kl, x_kt, Yl, Yt);
    end
    NRMSE(step,stepTrial)=N; NRMSE_C(step,stepTrial)=NC;
end

NRMSE = [mean(NRMSE,2,'omitnan') std(NRMSE,0,2,'omitnan') trial-sum(isnan(NRMSE),2) paramsSet NaN(searchParams,1) NRMSE];
NRMSE_C = [mean(NRMSE_C,2,'omitnan') std(NRMSE_C,0,2,'omitnan') trial-sum(isnan(NRMSE_C),2) paramsSet  NaN(searchParams,1) NRMSE_C];

[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,1));
bestNRMSE_C_ofBestNRMSE = NRMSE_C(bestNRMSEIndex,1);

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'TDRC=', Model, '_', Task, num2str(order), '_biasCheck=', num2str(biasCheck),'_inputCheck=', num2str(inputCheck), ...
    '_c=', num2str(cMin),'-', num2str(cMax), '_p=', num2str(pMin),'-', num2str(pMax),...
    '_gamma=', num2str(gammaInit),'-', num2str(gammaFin), '_trial=', num2str(trial), '.mat'), '-v7.3');
end