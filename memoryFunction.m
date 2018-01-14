function  [NRMSE_pinv, l_start, step_half, N,  x_kl, x_kt, ul, ut, data, l_interval] = memoryFunction(lin, eta, c1, gamma)
% 非線形時間遅れによるRC

%% 初期設定
% data*step_allの行列のデータを入力
data = 1; % 1ステップあたりの入力データ数
step_half = 3200; % 初期値は3200
step_all = 2*step_half; % データ長=ステップ数
theta = 0.01; % [s] 0.2 0.01

tau = 80; %[s] 400*theta
N = tau/theta; % delayのstep数
l_num = 400;
% lin = 1;
% eta = 0.5; % 1.8 2 1.3
% c1 =-0.1; %0=線形モデル、-0.1=非線形モデル
% gamma =0.5;
taskDelay = 30;

%データの長さ
step_half =step_half*N;
step_all = step_all*N;




%% 目標データの生成
data_length = step_half/N+200; % 少し長めに準備
no_inf=1;

%格納変数・初期値
Yt_l = zeros(data,data_length);
Yt_t = zeros(data,data_length);
while no_inf>0
%     ul = 0.5*rand(data,data_length); % 入力データを学習とテストで作り変える
%     ut = 0.5*rand(data,data_length);
%     ul = [ones(1,100),zeros(1,data_length-100)];
    ul = normrnd(0,1,[data,data_length]); % 入力データを学習とテストで作り変える
    ut = normrnd(0,1,[data,data_length]);
    
    %     save('ul', 'ul');
    %     save('ut', 'ut');
    %     load('ul.mat');
    %     load('ut.mat');
    
        for y_step=taskDelay+1:data_length
            Yt_l(:,y_step) = ul(:,y_step-taskDelay);
            Yt_t(:,y_step) = ut(:,y_step-taskDelay);
        end
    
    no_inf=sum(sum(isnan(Yt_l)))+sum(sum(isinf(Yt_l)))+sum(sum(isnan(Yt_t)))+sum(sum(isinf(Yt_t)));
end


%% 周期の算出
X_step = step_all*2; %数が足りる最小の長さで
% 周期固定の場合
i_period = N;
index=step_all/N;
step = (1:step_all/N+1)';
p_start = (step-1)*i_period+1;
p_interval(:,:) = p_start(2:index,:) - p_start(1:index-1,:);


%% Reservoirへのデータの挿入
p_begin = 1; % 15
W = 1; % 入力の波の数 p_intervalでの波長の合計が400を超えるようにとる
p_max = max(p_interval(p_begin:end,1));
% load('preM.mat');
preM = -0.1+0.2*round(rand(data,p_max*W*1)); % Mask
% save('preM','preM');

% 長めの入力を作成
Jl = zeros(length(preM),data_length, data); %行がマスク，列が入力uの2次元配列
Jt = zeros(length(preM),data_length, data);
for step = 1:data_length
    Jl(:,step) = ul(:,step).*preM;
    Jt(:,step) = ut(:,step).*preM;
end

% リザーバのデータを毎ステップとって，周期の始点を検出しながらMGの計算
Xl = zeros(1,X_step/2);
Xt = zeros(1,X_step/2);
X_index_l = ones(data_length,2,data); % 入力した長さを記録
X_index_t = ones(data_length,2,data);

X_tau = 0; % Xl(:,1)
X_tau2 =0; %Xt(:,1)
u_step_l = 1;
u_step_t = 1;
M_step_l = 1;
M_step_t = 1;
for step = 2:N
    M_step_l = M_step_l +1;
    dxl = dif(Xl(:,step-1),X_tau, lin, eta, c1, gamma, Jl(M_step_l,u_step_l));
    Xl(:,step) = Xl(:,step-1) + dxl*theta;
    
    M_step_t = M_step_t +1;
    dxt = dif(Xt(:,step-1),X_tau2, lin, eta, c1, gamma, Jt(M_step_t,u_step_t));
    Xt(:,step) = Xt(:,step-1) + dxt*theta;
end

while u_step_l < data_length
    if step >= p_start(u_step_l+1,1)-1
        X_index_l(u_step_l,2) = M_step_l;
        X_index_l(u_step_l+1, 1) = X_index_l(u_step_l, 1)+M_step_l;
        M_step_l = 0;
        u_step_l =u_step_l+1;
    end
    step = step+1;
    M_step_l = M_step_l +1;
    X_tau = Xl(:,step-N);
    dxl = dif(Xl(:,step-1), X_tau, lin, eta, c1, gamma, Jl(M_step_l,u_step_l));
    Xl(:,step) = Xl(:,step-1) + dxl*theta;
end

step = N;
while u_step_t < data_length
    if step >= p_start(u_step_t+1,1)-1
        X_index_t(u_step_t,2) = M_step_t;
        X_index_t(u_step_t+1, 1) = X_index_t(u_step_t, 1)+M_step_t;
        M_step_t = 0;
        u_step_t =u_step_t+1;
    end
    step = step+1;
    M_step_t = M_step_t +1;
    X_tau2 = Xt(:,step-N);
    dxt = dif(Xt(:,step-1), X_tau2, lin, eta, c1, gamma, Jt(M_step_t,u_step_t));
    Xt(:,step) = Xt(:,step-1) + dxt*theta;
end

%% 学習データから重みを求める memoryFunction
l_interval = N;
l_start = 100 + 50; % X(X0)が安定したところから学習スタート
t_start = l_start; % テストスタート
x_kl = zeros(l_num, data_length);
x_kt = zeros(l_num, data_length);
for k_step = 1:data_length-1
    x_kl(:,k_step) = Xl(:,X_index_l(k_step,1):l_interval/l_num:X_index_l(k_step,1)+l_interval-1);
    x_kt(:,k_step) = Xt(:,X_index_t(k_step,1):l_interval/l_num:X_index_t(k_step,1)+l_interval-1);
end


C = zeros(l_num,l_num);
p = zeros(l_num,data);
for step = l_start:l_start+step_half/N-1
    C = C + x_kl(:,step)*x_kl(:,step)';
    p = p + Yt_l(:,step)*x_kl(:,step);
end
C = C/(step_half/N);
p = p/(step_half/N);

Wout = zeros(data,l_num);
Wout(1,:) = pinv(C)*p;
%
% % Y_pinv_C = zeros(step_half/N,3);
% % Y_pinv_C(:,1:2) = [(Yt_l(:,l_start:l_start+step_half/N-1))' (Wout(1,:)*x_kl(:,l_start:l_start+step_half/N-1))'];
% % Y_pinv_C(:,3) =  Y_pinv_C(:,1).*Y_pinv_C(:,2);
% % mf_C = mean(Y_pinv_C(:,3));
%
Y_pinv_C = zeros(step_half/N,3);
Y_pinv_C(:,1:2) = [(Yt_l(:,l_start:l_start+step_half/N-1))' (Wout(1,:)*x_kl(:,l_start:l_start+step_half/N-1))'];
Y_pinv_C(:,3) =  (Y_pinv_C(1:end,1)-Y_pinv_C(1:end,2)).^2;
sigma_ytC=std(Y_pinv_C(1:end,1));
NRMSE_pinv_C =((1/(step_half/N))*sum(Y_pinv_C(1:end,3))/(sigma_ytC^2))^0.5;
%
t_begin = t_start;
t_fin = t_start+step_half/N-1;

Y_pinv = zeros(step_half/N,3);
Y_pinv(:,1:2) = [(Yt_t(:, t_begin:t_fin))' (Wout(1,:)*x_kt(:, t_begin:t_fin))'];
% Y_pinv(:,3) =  Y_pinv(:,1).*Y_pinv(:,2);
% mf = mean(Y_pinv(:,3));
Y_pinv(:,3) =  (Y_pinv(1:end,1)-Y_pinv(1:end,2)).^2;
sigma_yt=std(Y_pinv(1:end,1));
NRMSE_pinv =((1/ (step_half/N))*sum(Y_pinv(1:end,3))/(sigma_yt^2))^0.5;
%
% %% 学習データから重みを求める NRMSE
% %学習に必要な数だけまとめる
% % l_start = 50; % X(X0)が安定したところから学習スタート
% % t_start = l_start; % テストスタート
% % l_interval = N;
% % x_kl = zeros(l_interval, data_length);
% % x_kt = zeros(l_interval, data_length);
% % for k_step = 1:data_length-1
% %     x_kl(:,k_step) = Xl(:,X_index_l(k_step,1):X_index_l(k_step,1)+l_interval-1);
% %     x_kt(:,k_step) = Xt(:,X_index_t(k_step,1):X_index_t(k_step,1)+l_interval-1);
% % end
%
% % %% 疑似逆行列テスト
% % % 疑似逆行列
% % Wout = zeros(data,l_interval);
% % Wout(1,:) =Yt_l(1,l_start:l_start+step_half/N-1)*pinv(x_kl(:,l_start:l_start+step_half/N-1));
% %
% % t_begin = t_start;
% % t_fin = t_start+step_half/N-1;
% %
% % Y_pinv_C = zeros(step_half/N,3);
% % Y_pinv_C(:,1:2) = [(Yt_l(:,l_start:l_start+step_half/N-1))' (Wout(1,:)*x_kl(:,l_start:l_start+step_half/N-1))'];
% % Y_pinv_C(:,3) =  (Y_pinv_C(1:end,1)-Y_pinv_C(1:end,2)).^2;
% % sigma_ytC=std(Y_pinv_C(1:end,1));
% % NRMSE_pinv_C =((1/(step_half/N))*sum(Y_pinv_C(1:end,3))/(sigma_ytC^2))^0.5;
% %
% % Y_pinv = zeros(step_half/N,3);
% % Y_pinv(:,1:2) = [(Yt_t(:, t_begin:t_fin))' (Wout(1,:)*x_kt(:, t_begin:t_fin))'];
% % Y_pinv(:,3) =  (Y_pinv(1:end,1)-Y_pinv(1:end,2)).^2;
% % sigma_yt=std(Y_pinv(1:end,1));
% % NRMSE_pinv =((1/ (step_half/N))*sum(Y_pinv(1:end,3))/(sigma_yt^2))^0.5;
% %
% %
% % %% 特異値分解テスト
% % % X_tmp=x_kl(:,l_start:l_start + step_half/N-1);
% % % [U,S,V] = svd(X_tmp);
% % % [A,B]=size(S);
% % % for i=1:A
% % %     for j=1:B
% % %         if abs(S(i,j)) >1*10^-17 % ここを変えたらいいかも
% % %             S(i,j)=1/S(i,j);
% % %             %S(i,j)=S(i,j)/( S(i,j)^2+reg^2); %Tikhonov Regularization
% % %         else
% % %             S(i,j)=0;
% % %         end
% % %     end
% % % end
% % % Wout2 =Yt_l(1,l_start:l_start+step_half/N-1)*(V*S'*U');
% % % % load('Wout2.mat');
% % %
% % % Y_SVD_C = zeros(step_half/N,3);
% % % Y_SVD_C(:,1:2) = [(Yt_l(:,t_start:t_start+step_half/N-1))' (Wout2(1,:)*x_kl(:,t_start:t_start+step_half/N-1))'];
% % % Y_SVD_C(:,3) =  (Y_SVD_C(1:end,1)-Y_SVD_C(1:end,2)).^2;
% % % sigma_ytC=std(Y_SVD_C(1:end,1));
% % % NRMSE_SVD_C =((1/(step_half/N))*sum(Y_SVD_C(1:end,3))/(sigma_ytC^2))^0.5;
% % %
% % % Y_SVD = zeros(step_half/N,3);
% % % Y_SVD(:,1:2) = [(Yt_t(:, t_begin:t_fin))' (Wout2(1,:)*x_kt(:, t_begin:t_fin))'];
% % % Y_SVD(:,3) =  (Y_SVD(1:end,1)-Y_SVD(1:end,2)).^2;
% % % sigma_yt2=std(Y_SVD(1:end,1));
% % % NRMSE_SVD =((1/ (step_half/N))*sum(Y_SVD(1:end,3))/(sigma_yt2^2))^0.5;
% % % Xl;
end

%% モデル定義
function dx = dif(X, X_tau, lin, eta, c1, gamma, J)
dx = -lin*X + eta*X_tau + c1*X^3 + gamma*J;
end
