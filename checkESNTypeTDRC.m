function [bestNRMSE] = checkESNTypeTDRC(ang)
load('\\NAS08C093\Public\kinoshita\test_nas\20180726MaxEigen=2_check\201807271838TDRC=NL_NARMA10_biasCheck=0_inputCheck=0_c=-2.5-1.5_p=3-3_gamma=-2-1_trial=100.mat')
[best,b]=min(NRMSE(:,6));

a = NRMSE(b,1); b1=NRMSE(b,2); b2=b1/2; gamma = NRMSE(b,5); oneDim=200;
% ang = pi/4;
P = [cos(ang) -sin(ang); sin(ang) cos(ang)];
eigMat = [b1 0; 0 b2];
W = inv(P)*eigMat*P;

[bestNRMSE] = ESNTypeTDRC(100,1,theta, oneDim, biasCheck, inputCheck, a, gamma, W,10,0,ang);
end