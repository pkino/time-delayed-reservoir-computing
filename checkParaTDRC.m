clear

load('\\NAS08C093\Public\kinoshita\allData\simulationData\20180726MaxEigen=2_check\201807271838TDRC=NL_NARMA10_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-2-1_trial=100.mat')
[a,b]=min(NRMSE(:,6))
paraParam(1,:) = [theta, 400, biasCheck, inputCheck, NRMSE(b,1:5)]

load('\\NAS08C093\Public\kinoshita\allData\simulationData\20180726MaxEigen=2_check\201901281342TDRC=NL_NARMA30_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-2-1_trial=100.mat')
[a,b]=min(NRMSE(:,6))
paraParam(2,:) = [theta, 400, biasCheck, inputCheck, NRMSE(b,1:5)]

paraTDRC(100,2,paraParam,10,0);
paraTDRC(100,2,paraParam,30,0);