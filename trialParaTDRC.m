function trialParaTDRC(order,paraNum)
dbstop if error
trial=100; Task='NARMA'; theta=0.2; learnDimension=400/paraNum; biasCheck=0; inputCheck=0; p=3;
diffMask=1;

cNum = 3;
cMin = -2.5;
cMax = 1.5;
cData = logspace(cMin,cMax,cNum);

gammaNum = 3;
gammaMin = -2;
gammaMax = 1;
gammaData = logspace(gammaMin, gammaMax, gammaNum);

% マップの固有値の範囲を指定
eigMin=-1.5; eigMax=-6; eigNum=5; gapMax=-2; gapMin=-5.5; gapNum=5; % eigMin=-1.5; eigMax=-9; gapMax=-2; gapMin=-8;
[a_mat, b_mat] = find_ab(eigMin, eigMax, eigNum, gapMax, gapMin, gapNum);

RCLen=1500; % 最終の学習の長さ
dataLen = RCLen + 500;
dataGen = str2func(strcat('dataGenerator_',Task));

singleSearchNum = 4;
singleSearchParams = gapNum*eigNum*cNum*gammaNum;
singleParamsSet = [];
singleParamsSetIndex = (1:singleSearchParams)';

for step_gap = 1:gapNum
    for step_eig = 1:eigNum
        for step_c = 1:cNum
            for step_g = 1:gammaNum
                singleParamsSet = [singleParamsSet; a_mat(step_gap,step_eig), b_mat(step_gap,step_eig), ...
                    cData(step_c), gammaData(step_g)];
            end
        end
    end
end

paramsSetIndex = singleParamsSetIndex;
for step = 1:paraNum-1
    paramsSetIndex = [repelem(singleParamsSetIndex, length(paramsSetIndex),1) repmat(paramsSetIndex, singleSearchParams,1)];
end
for step = 1:paraNum-1
    paramsSetIndex= paramsSetIndex(paramsSetIndex(:,step)<=paramsSetIndex(:,step+1),:);
end
paramsSet = [];
searchParams = length(paramsSetIndex);
parfor step =1:searchParams
    tempParamsSet = [];
    for step2 = 1:paraNum
        tempParamsSet = [tempParamsSet singleParamsSet(paramsSetIndex(step,step2),:)];
    end
    paramsSet = [paramsSet; tempParamsSet];
end

NRMSE = zeros(searchParams,trial); NRMSE_C = zeros(searchParams,trial);
parfor step = 1:searchParams
    for stepTrial = 1:trial
        %% 入出力データ作成
        [ul, ~, ut, Yl, ~, Yt] = dataGen(dataLen,order);
        
        
        %% リザーバ計算
        para_x_kl = cell(paraNum); para_x_kt = cell(paraNum);
        for stepPara = 1:paraNum
            [para_x_kl{stepPara}, para_x_kt{stepPara}] = ...
                timeDelayReservoir(stepTrial, ul, ut, theta, learnDimension, ...
                biasCheck, inputCheck, paramsSet(step,(stepPara-1)*singleSearchNum+1), ...
                paramsSet(step,(stepPara-1)*singleSearchNum+2),  ...
                paramsSet(step,(stepPara-1)*singleSearchNum+3), p, ...
                paramsSet(step,(stepPara-1)*singleSearchNum+4));
        end
        
        %% 学習とテスト
        try
            [NRMSE(step,stepTrial), NRMSE_C(step,stepTrial)] = TDRC(RCLen, vertcat(para_x_kl{:}), vertcat(para_x_kt{:}), Yl, Yt);
        catch
            NRMSE(step,stepTrial) =  NaN; NRMSE_C(step,stepTrial) =  NaN;
        end
    end
end

NRMSE = [mean(NRMSE,2,'omitnan') std(NRMSE,0,2,'omitnan') trial-sum(isnan(NRMSE),2) paramsSet NaN(searchParams,1) NRMSE];
NRMSE_C = [mean(NRMSE_C,2,'omitnan') std(NRMSE_C,0,2,'omitnan') trial-sum(isnan(NRMSE_C),2) paramsSet  NaN(searchParams,1) NRMSE_C];

[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,1));
bestNRMSE_C_ofBestNRMSE = NRMSE_C(bestNRMSEIndex,1);

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'paraTDRC=',num2str(paraNum), '_', Task ,num2str(order), '_c=', num2str(cMin), '-', ...
    num2str(cMax), '_gamma=', num2str(gammaMin), '-', num2str(gammaMax), '_diffMask=', num2str(diffMask) ,'.mat'), '-v7.3');
end
