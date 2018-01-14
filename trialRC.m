function trialRC(taskDelay, trial, initNum)
load(strcat('Res_',num2str(initNum),'.mat'),'eig_num','gap_num','g_num');
dataNum = eig_num*gap_num*g_num;
REC_NARMA = NaN(dataNum+1,1+2+1+trial);
data_NARMA = NaN(g_num, gap_num, eig_num);


for step = 1:trial
    load(strcat('Res_',num2str(initNum+step-1),'.mat'));
    gap_num=gap_num;
    g_num=g_num;
    
    parfor step_eig =1:eig_num
        for step_gap =1:gap_num
            for step_g = 1:g_num
                [step_half, N, l_num,x_kl, x_kt, ul, ut, l_start] = Res_data{:, step_g, step_gap, step_eig};
                data_NARMA(step_g,step_gap,step_eig) = RC_MF(taskDelay, step_half, N, l_num,x_kl, x_kt, ul, ut, l_start);
            end
        end
    end
    REC_NARMA(1:end-1,4+step) = reshape(data_NARMA,[],1);
end

if g_num==1
    g_data = g_init;
else
    g_data = logspace(g_init,g_fin,g_num)';
end
REC_NARMA(1:end-1,1) = repmat(g_data,eig_num*gap_num,1);
REC_NARMA(1:end-1,2) = mean(REC_NARMA(1:end-1,5:end),2);
REC_NARMA(1:end-1,3) = std(REC_NARMA(1:end-1,5:end),0,2);

Date = datestr(datetime('now'),'yyyymmddHHMM');
save(strcat(Date,'uniform_', Model,'_MF',num2str(taskDelay),'_NARMA','fixedEigenDistance_', 'gamma_expPortion=', num2str(g_init),'-', num2str(g_fin), '.mat'), 'REC_NARMA','eigMin', 'eigMax', 'gapMax', 'gapMin','g_init', 'g_fin','eig_num','gap_num','g_num','taskDelay');
end
