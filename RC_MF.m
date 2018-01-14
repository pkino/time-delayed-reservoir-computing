function [NRMSE_pinv] = RC_MF(taskDelay, step_half, N, l_num,x_kl, x_kt, ul, ut, l_start)
Yt_l = [zeros(1, taskDelay), ul];
Yt_t = [zeros(1, taskDelay), ut];

C = zeros(l_num,l_num);
p = zeros(l_num,1);
for l_step = l_start:l_start+step_half/N-1
    C = C + x_kl(:,l_step)*x_kl(:,l_step)';
    p = p + Yt_l(:,l_step)*x_kl(:,l_step);
end
C = C/(step_half/N);
p = p/(step_half/N);


Wout = zeros(1,l_num);
Wout(1,:) = pinv(C)*p;

t_begin = l_start;
t_fin = l_start+step_half/N-1;

% Y_pinv_C = zeros(step_half/N,3);
% Y_pinv_C(:,1:2) = [(Yt_l(:,l_start:l_start+step_half/N-1))' (Wout(1,:)*x_kl(:,l_start:l_start+step_half/N-1))'];
% Y_pinv_C(:,3) =  (Y_pinv_C(1:end,1)-Y_pinv_C(1:end,2)).^2;
% sigma_ytC=std(Y_pinv_C(1:end,1));
% NRMSE_pinv_C =((1/(step_half/N))*sum(Y_pinv_C(1:end,3))/(sigma_ytC^2))^0.5;

Y_pinv = zeros(step_half/N,3);
Y_pinv(:,1:2) = [(Yt_t(:, t_begin:t_fin))' (Wout(1,:)*x_kt(:, t_begin:t_fin))'];
Y_pinv(:,3) =  (Y_pinv(1:end,1)-Y_pinv(1:end,2)).^2;
sigma_yt=std(Y_pinv(1:end,1));
NRMSE_pinv =((1/ (step_half/N))*sum(Y_pinv(1:end,3))/(sigma_yt^2))^0.5;

end

