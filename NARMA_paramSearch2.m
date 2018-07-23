function NARMA_paramSearch2(order, fixedC)

load('\\WORKSTATION2013\test\20180628smallMap_p=3_Linear_NonLinearModel\201806280847TDRC=NL_NARMA10_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-4-1.5_trial=20.mat')
paramMap1 = paramSearchFromMap(NRMSE, eigNum, gapNum, gammaNum, cNum, fixedC);

load('\\WORKSTATION2013\test\20180628smallMap_p=3_Linear_NonLinearModel\201806280947TDRC=NL_NARMA30_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-4-1.5_trial=20.mat')
paramMap2 = paramSearchFromMap(NRMSE, eigNum, gapNum, gammaNum, cNum, fixedC);


%　各変数定義しなおし
trial=1; initNo=1; gammaCheck=0;

NRMSE_num = (eigNum*gapNum)^2;
NRMSE = NaN(NRMSE_num, 4+3+trial); NRMSE_C = NaN(NRMSE_num, 4+3+trial); seedREC = cell(NRMSE_num,2);

recIndex=0;
for map1Gap = 1:gapNum
    for map1Eig = 1:eigNum
        for map2Gap = 1:gapNum
            for map2Eig = 1:eigNum
                paraParam = [theta, 400, biasCheck, inputCheck,paramMap1{map1Gap, map1Eig}; ...
                    theta, 400, biasCheck, inputCheck, paramMap2{map2Gap, map2Eig}];
                
                recIndex = recIndex+1;
                NRMSE(recIndex,1:4) = [map1Gap,map1Eig,map2Gap,map2Eig];
                NRMSE_C(recIndex,1:4) = [map1Gap,map1Eig,map2Gap,map2Eig];
                seedREC{recIndex,1} = [map1Gap,map1Eig,map2Gap,map2Eig];
                [NRMSE(recIndex,5:end), NRMSE_C(recIndex,5:end), seedREC{recIndex,2}] = trialParaTDRC(trial,initNo,2,paraParam,order,gammaCheck,1);
            end
        end
    end
end

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'searchParamSet_paraTDRC',num2str(2), '_NARMA',num2str(order), '_fiexdC=', num2str(fixedC),'_gammaCheck=', num2str(gammaCheck), '_trial=', num2str(trial), '.mat'), '-v7.3');
end

function [paramMap] = paramSearchFromMap(NRMSE, eigNum, gapNum, gammaNum, cNum, fixedC)
data = NRMSE;
saveNum = 5;

cIndex = 3;
[bestPerform, bestIndex] = min(data(:,saveNum+1));
bestC = data(bestIndex, cIndex);

horizontal = eigNum;
vertical = gapNum;
other = gammaNum*cNum;

paramMap = cell(vertical, horizontal);
for step = 1:vertical
    for step2 = 1:horizontal
        param = data((step-1)*vertical*other+(step2-1)*other+1:(step-1)*vertical*other+(step2-1)*other+other, 1:saveNum+2);
        
        if fixedC == 1
            % c固定
            param = param(param(:,cIndex) == bestC,:);
        end
        
        
        [gridMin, miniIndex]  = min(param(:,end-1));
        paramMap{step,step2} = param(miniIndex,1:saveNum);
        
        % 各値が知りたいときはコメントアウトを外す
        %         miniIndex
        %         c=plotter(miniIndex,3)
    end
end

end
