clear

load('\\WORKSTATION2013\test\20180628smallMap_p=3_Linear_NonLinearModel\201806280847TDRC=NL_NARMA10_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-4-1.5_trial=20.mat')
[a,b]=min(NRMSE(:,6))
paraParam(1,:) = [theta, 400, biasCheck, inputCheck, NRMSE(b,1:5)]

load('\\WORKSTATION2013\test\20180628smallMap_p=3_Linear_NonLinearModel\201806280947TDRC=NL_NARMA30_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-4-1.5_trial=20.mat')
[a,b]=min(NRMSE(:,6))
paraParam(2,:) = [theta, 400, biasCheck, inputCheck, NRMSE(b,1:5)]

trialParaTDRC(100,1,1,paraParam,10,0,0);
trialParaTDRC(100,1,2,paraParam,30,0,0);