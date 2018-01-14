function trialRes
trial = 3;
initNum = 1;
Model = 'L';
if Model == 'L'
    c1 = 0;
    g_init = 3;
    g_num = 3;
    g_fin = -3;
    g_data = logspace(g_init, g_fin, g_num)';
elseif Model == 'NL'
    c1= -0.1;
    g_init = 0.1;
    g_fin = 0.1;
    g_inter = 0.1;
    g_num = int32((g_fin-g_init)/g_inter+1);
    g_data = (g_init:g_inter:g_fin)';
else
    error('ƒ‚ƒfƒ‹‚ğ³‚µ‚­‘I‘ğ‚µ‚Ä‚­‚¾‚³‚¢');
end

eigMin=-1.5; eigMax=-9; gapMax=-2; gapMin=-8;
[a_mat, b_mat] = find_ab(eigMin, eigMax, gapMax, gapMin);

gap_num = size(a_mat,1);
eig_num = size(a_mat,2);

Date = datestr(datetime('now'),'yyyymmddHHMM');
for step = 1:trial
    Res_data = cell(8, g_num, eig_num, gap_num);
    rngNum = initNum+step-1
    parfor step_eig =1:eig_num
        for step_gap =1:gap_num
            a = a_mat(step_gap,step_eig);
            b = b_mat(step_gap,step_eig);
            for step_g = 1:g_num
                gamma = g_data(step_g,1);
                rng(rngNum);
                [Res_data{:,step_g,step_gap,step_eig}] = Res(a,b,c1,gamma);
            end
        end
    end
    save(strcat('Res_', num2str(rngNum)),'Res_data','gap_num', 'eig_num', 'g_num', 'eigMin', 'eigMax', 'gapMax', 'gapMin','g_init', 'g_fin','Date','Model','-v7.3');
end
end


