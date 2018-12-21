function trialSearchSeriTDRC(order,seriNum)
dbstop if error
trial=2; initNo=1; Task='NARMA'; theta=0.2; learnDimension=400; biasCheck=0; inputCheck=0; p=3;
% seriNum=2; 
diffMask=1;

cNum = 3;
cMin = -2.5;
cMax = 1.5;
cData = logspace(cMin,cMax,cNum);

gammaNum = 3;
gammaMin = -2;
gammaMax = 1;
gammaData = logspace(gammaMin, gammaMax, gammaNum);

eigMin=-1.5; eigMax=-6; eigNum=5; gapMax=-2; gapMin=-5.5; gapNum=5; % eigMin=-1.5; eigMax=-9; gapMax=-2; gapMin=-8;
[a_mat, b_mat] = find_ab(eigMin, eigMax, eigNum, gapMax, gapMin, gapNum);

RCLen=1500; % 最終の学習の長さ
oneSeriSteps = 200;  dataLen = RCLen+seriNum*oneSeriSteps; %oneSteps=l_Start(TDRC)と設定
seed_status = 1; % seed_dataGen - seed_mask
dataGen = str2func(strcat('dataGenerator_',Task));

saveNum = 2;
singleSearchNum = 4;
searchNum = singleSearchNum*seriNum;
singleSearchParams = gapNum*eigNum*cNum*gammaNum;
searchParams = singleSearchParams^seriNum;
singleParamsSet = NaN(singleSearchParams,singleSearchNum);

index = 1;
for step_gap = 1:gapNum
    for step_eig = 1:eigNum
        for step_c = 1:cNum
            for step_g = 1:gammaNum
                singleParamsSet(index, :) = [a_mat(step_gap,step_eig), b_mat(step_gap,step_eig), ...
                    cData(step_c), gammaData(step_g)];
                index = index + 1;
            end
        end
    end
end

paramsSet = singleParamsSet;
for step = 1:seriNum-1
    paramsSet = [repelem(singleParamsSet, length(paramsSet),1) repmat(paramsSet, singleSearchParams,1)];
end

NRMSE = zeros(searchParams,trial); NRMSE_C = zeros(searchParams,trial);
parfor step = 1:searchParams
    for stepTrial = 1:trial
        %% 入出力データ作成
        seed_dataGen = stepTrial+initNo-1;
        [ul, ~, ut, Yl, ~, Yt] = dataGen(dataLen,seed_dataGen,order);
        
        %% 1つ目のRC
        cycLen = dataLen;
        seed_no = seed_dataGen + seed_status;
        [x_kl, x_kt] ...
            = timeDelayReservoir(seed_no, ul, ut, theta, learnDimension, biasCheck, inputCheck, ...
            paramsSet(step,1), paramsSet(step,2), paramsSet(step,3), p, paramsSet(step,4));
        cycLen = cycLen-oneSeriSteps;
        [N,NC,Ylp,Ytp] = TDRC(cycLen, x_kl, x_kt, Yl, Yt);
        
        %% 2つ目以降のRC
        try % 1つ目のRCでYlp, YtpがNaNになる可能性があるため
            for stepSeri = 2:seriNum
                if diffMask==1
                    seed_mask = seed_no + stepSeri -1;
                else
                    seed_mask = seed_no;
                end
                
                %% リザーバ計算
                [x_kl, x_kt] ...
                    = timeDelayReservoir(seed_mask, Ylp', Ytp', theta, learnDimension, biasCheck, inputCheck, ...
                    paramsSet(step,(stepSeri-1)*singleSearchNum+1), paramsSet(step,(stepSeri-1)*singleSearchNum+2), ...
                    paramsSet(step,(stepSeri-1)*singleSearchNum+3), p, paramsSet(step,(stepSeri-1)*singleSearchNum+4));
                
                %% 学習とテスト
                Yl = Yl(oneSeriSteps:end-1); Yt=Yt(oneSeriSteps:end-1);  % 一つ前のTDRCのx_kl, x_ktのl_start=200のため
                cycLen = cycLen-oneSeriSteps;
                [N,NC,Ylp,Ytp] = TDRC(cycLen, x_kl, x_kt, Yl,Yt);
            end
            NRMSE(step,stepTrial)=N; NRMSE_C(step,stepTrial)=NC;
        catch
            NRMSE(step,stepTrial)=NaN; NRMSE_C(step,stepTrial)=NaN;
        end
    end
end

NRMSE = [mean(NRMSE,2,'omitnan') std(NRMSE,0,2,'omitnan') trial-sum(isnan(NRMSE),2) paramsSet NaN(searchParams,1) NRMSE];
NRMSE_C = [mean(NRMSE_C,2,'omitnan') std(NRMSE_C,0,2,'omitnan') trial-sum(isnan(NRMSE_C),2) paramsSet  NaN(searchParams,1) NRMSE_C];

[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,1));
bestNRMSE_C_ofBestNRMSE = NRMSE_C(bestNRMSEIndex,1);

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'seriTDRC=',num2str(seriNum), '_', Task ,num2str(order), '_c=', num2str(cMin), '-', ...
    num2str(cMax), '_gamma=', num2str(gammaMin), '-', num2str(gammaMax), '_diffMask=', num2str(diffMask) ,'.mat'), '-v7.3');
end