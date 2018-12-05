function trialSeriTDRC(seriNum, order, diffMask)
trial=100;

NRMSE = NaN(1, trial);
NRMSE_C = NaN(1, trial);

oneSteps = 500; dataLen = 1500+seriNum*oneSteps;
for stepTrial = 1:trial
    seed_dataGen = stepTrial;
    [ul, ~, ut, Yl, ~, Yt] = dataGenerator_modifiedNARMA(dataLen,seed_dataGen,order);
    
    
    cycLen = dataLen;
    seedStatus=1;
    seed_no = seed_dataGen + seedStatus;
    %% リザーバ計算
    [x_kl, x_kt] ...
        = timeDelayReservoir(seed_no, ul, ut, 0.2, 400, 0, 0, ...
        0.4, 0.3, 0.003, 3, 0.1);
    
    %% 学習とテスト
    cycLen = cycLen-500;
    [N,NC,Ylp,Ytp] = TDRC(cycLen, x_kl, x_kt, Yl, Yt);
    
    
    for step = 2:seriNum
        if diffMask==1
            seed_mask = seed_no + step -1;
        else
            seed_mask = seed_no;
        end
        
        %% リザーバ計算
        [x_kl, x_kt] ...
            = timeDelayReservoir(seed_mask, Ylp', Ytp', 0.2, 400, 0, 0, ...
            0.4, 0.3, 0.003, 3, 0.01);
        
        %% 学習とテスト
        Yl = Yl(200:end); Yt=Yt(200:end);
        cycLen = cycLen-oneSteps;
        [N,NC,Ylp,Ytp] = TDRC(cycLen, x_kl, x_kt, Yl,Yt);
    end
    NRMSE(1,stepTrial)=N; NRMSE_C(1,stepTrial)=NC;
end
NRMSE = [mean(NRMSE,'omitnan') std(NRMSE,0,2,'omitnan') trial-sum(isnan(NRMSE),2) NaN NRMSE];
NRMSE_C = [mean(NRMSE_C,'omitnan') std(NRMSE_C,0,2,'omitnan') trial-sum(isnan(NRMSE_C),2)  NaN NRMSE_C];

[bestNRMSE, bestNRMSEIndex] = min(NRMSE(:,1));
bestNRMSE_C_ofBestNRMSE = NRMSE_C(bestNRMSEIndex,1);
[bestNRMSE_C, bestNRMSE_CIndex] = min(NRMSE_C(:,1));

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'seriTDRC=',num2str(seriNum), '_modifiedNARMA',num2str(order), '_diffMask=', num2str(diffMask), '_trial=', num2str(trial), '.mat'), '-v7.3');
end