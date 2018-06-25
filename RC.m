function [NRMSE_pinv, NRMSE_pinv_C] = RC(step_half, N, l_num,x_kl, x_kt, Yl, Yt, l_start)
Wout = zeros(1,l_num);
Wout(1,:) = Yl(1,l_start:l_start+step_half/N-1)*pinv(x_kl(:,l_start:l_start+step_half/N-1));

t_begin = l_start;
t_fin = l_start+step_half/N-1;

Y_pinv_C = zeros(step_half/N,3);
Y_pinv_C(:,1:2) = [(Yl(:,l_start:l_start+step_half/N-1))' (Wout(1,:)*x_kl(:,l_start:l_start+step_half/N-1))'];
Y_pinv_C(:,3) =  (Y_pinv_C(1:end,1)-Y_pinv_C(1:end,2)).^2;
sigma_ytC=std(Y_pinv_C(1:end,1));
NRMSE_pinv_C =((1/(step_half/N))*sum(Y_pinv_C(1:end,3))/(sigma_ytC^2))^0.5;

Y_pinv = zeros(step_half/N,3);
Y_pinv(:,1:2) = [(Yt(:, t_begin:t_fin))' (Wout(1,:)*x_kt(:, t_begin:t_fin))'];
Y_pinv(:,3) =  (Y_pinv(1:end,1)-Y_pinv(1:end,2)).^2;
sigma_yt=std(Y_pinv(1:end,1));
NRMSE_pinv = ((1/ (step_half/N))*sum(Y_pinv(1:end,3))/(sigma_yt^2))^0.5;
end

