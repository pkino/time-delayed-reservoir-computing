function trialTDRC(trial, initNo, Model, Task, order, theta, learnDimension, biasCheck, inputCheck, ...
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

pNum = 1;
% pMin = 3;
pMax = pMin;
pData = linspace(pMin,pMax,pNum);

searchNum = 5;
saveNum = 2;
saveData = cell(trial, gammaNum, pNum, cNum, eigNum, gapNum, saveNum);
paramLength = gapNum*eigNum*cNum*pNum*gammaNum;
paramsSet = NaN(paramLength,searchNum+3+1);

RCLen = 1500;
dataLen = RCLen + 500;
seed_status = 1; % seed_dataGen - seed_mask
dataGen = str2func(strcat('dataGenerator_',Task));

parfor step_gap = 1:gapNum
    for step_eig = 1:eigNum
        for step_c = 1:cNum
            for step_p = 1:pNum
                for step_g = 1:gammaNum
                    for stepTrial = 1:trial
                        %% 入力・目標データの生成
                        seed_dataGen = stepTrial+initNo-1;
                        [ul, ~, ut, Yl, ~, Yt] = dataGen(dataLen,seed_dataGen,order);
                        seed_no = stepTrial+seed_status;
                        
                        %% リザーバ計算
                        [x_kl, x_kt] ...
                            = timeDelayReservoir(seed_no, ul, ut, theta, learnDimension, biasCheck, inputCheck, ...
                            a_mat(step_gap,step_eig), b_mat(step_gap,step_eig), cData(step_c), pData(step_p), gammaData(step_g));
                        
                        %% 学習とテスト
                        [saveData{stepTrial,step_g, step_p, step_c, step_eig, step_gap,:}] = TDRC(RCLen, x_kl, x_kt, Yl, Yt);
                    end
                end
            end
        end
    end
end
saveData = permute(reshape(saveData,trial,[],saveNum),[2 1 3]);

index = 1;
for step_gap = 1:gapNum
    for step_eig = 1:eigNum
        for step_c = 1:cNum
            for step_p = 1:pNum
                for step_g = 1:gammaNum
                    paramsSet(index, 1:searchNum) = [a_mat(step_gap,step_eig), b_mat(step_gap,step_eig), ...
                        cData(step_c), pData(step_p), gammaData(step_g)];
                    index = index + 1;
                end
            end
        end
    end
end

REC = zeros(paramLength, searchNum+3+1+trial,saveNum);
for stepSaveNum = 1:saveNum
    paramsSet(:,searchNum+1) = mean(cell2mat(saveData(:,:,stepSaveNum)),2,'omitnan');
    paramsSet(:,searchNum+2) = std(cell2mat(saveData(:,:,stepSaveNum)),0,2,'omitnan');
    paramsSet(:,searchNum+3) = trial - sum(isnan(cell2mat(saveData(:,:,stepSaveNum))),2);
    REC(:,:,stepSaveNum) =  horzcat(paramsSet,cell2mat(saveData(:,:,stepSaveNum)));
end

NRMSE = REC(:,:,1);
NRMSE_C = REC(:,:,2);

[bestNRMSE_C, bestNRMSE_CIndex] = min(NRMSE_C(:,searchNum+1));
[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,searchNum+1));

clear REC saveData

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'TDRC=', Model, '_NARMA',num2str(order), '_biasCheck=', num2str(biasCheck),'_inputCheck=', num2str(inputCheck), ...
    '_c=', num2str(cMin),'-', num2str(cMax), '_p=', num2str(pMin),'-', num2str(pMax),...
    '_gamma=', num2str(gammaInit),'-', num2str(gammaFin), '_trial=', num2str(trial), '.mat'), '-v7.3');
end