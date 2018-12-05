function [NRMSE_pinv, NRMSE_pinv_C, Ylp, Ytp] = TDRC(RCLen, x_kl, x_kt, Yl, Yt)
try
    %% 学習
    l_start = 200; % Xが安定したところから学習スタート
    
    Wout = Yl(1,l_start:l_start+RCLen-1)*pinv(x_kl(1:end,l_start:l_start+RCLen-1));
    
    %         Wout = (lasso( x_kl(:,l_start:l_start+RCLen-1)', Yl(1,l_start:l_start+RCLen-1)))';
    
    %     Wout =  Yl(1,l_start:l_start+RCLen-1)*x_kl(3:end,l_start:l_start+RCLen-1)'*inv(x_kl(3:end,l_start:l_start+RCLen-1)*x_kl(3:end,l_start:l_start+RCLen-1)');
    %     Wout = ridge(Yl(1,l_start:l_start+RCLen-1)',x_kl(1:end,l_start:l_start+RCLen-1)',1e-20,0)';
    
    t_begin = l_start;
    t_fin = l_start+RCLen-1;
    
    Y_pinv_C = zeros(RCLen,3);
    Ylp = (Wout(1,:)*x_kl(1:end,l_start:l_start+RCLen-1))';
    %     Ylp = (Wout(1,:)*[ones(1,RCLen);x_kl(1:end,l_start:l_start+RCLen-1)])';
    Y_pinv_C(:,1:2) = [(Yl(:,l_start:l_start+RCLen-1))' Ylp];
    Y_pinv_C(:,3) =  (Y_pinv_C(1:end,1)-Y_pinv_C(1:end,2)).^2;
    sigma_ytC=std(Y_pinv_C(1:end,1));
    NRMSE_pinv_C =((1/(RCLen))*sum(Y_pinv_C(1:end,3))/(sigma_ytC^2))^0.5;
    
    %% テスト
    Y_pinv = zeros(RCLen,3);
    Ytp = (Wout(1,:)*x_kt(1:end, t_begin:t_fin))';
    %     Ytp = (Wout(1,:)*[ones(1,RCLen);x_kt(1:end, t_begin:t_fin)])';
    Y_pinv(:,1:2) = [(Yt(:, t_begin:t_fin))' Ytp];
    Y_pinv(:,3) =  (Y_pinv(1:end,1)-Y_pinv(1:end,2)).^2;
    sigma_yt=std(Y_pinv(1:end,1));
    NRMSE_pinv = ((1/ (RCLen))*sum(Y_pinv(1:end,3))/(sigma_yt^2))^0.5;
catch
    NRMSE_pinv=NaN; NRMSE_pinv_C=NaN;Ylp=NaN; Ytp=NaN;
end
end

