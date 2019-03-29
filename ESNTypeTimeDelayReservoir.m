function  [x_kl, x_kt] = ...
    ESNTypeTimeDelayReservoir(seed_no, ul, ut, theta, learnDimension, ...
    biasCheck, inputCheck, a, gamma, W)
% step_half, N, learnDimension, x_kl, x_kt, ul, ut, l_start

%% 初期設定
% data*step_allの行列のデータを入力
data = 1; % 1ステップあたりの入力データ数
% theta = 0.01; % [s] 0.2 0.01

tau = 80; %[s] 400*theta
N = tau/theta; % delayのstep数
% lin = 1;
% eta = 0.5; % 1.8 2 1.3
% c1 =-0.1; %0=線形モデル、-0.1=非線形モデル
% gamma =0.5;

%データの長さ
dataLen = length(ul); % 少し長めに準備


%% 周期の算出
% 周期固定の場合
i_period = N;
index=dataLen;
step = (1:dataLen+1)';
p_start = (step-1)*i_period+1;
p_interval(:,:) = p_start(2:index,:) - p_start(1:index-1,:);


%% Reservoirへのデータの挿入
p_begin = 1; % 15
p_max = max(p_interval(p_begin:end,1));
seed_mask = seed_no;
rng(seed_mask,'twister');
maskCheck = 1;
while maskCheck > 0
    preM2 = -0.1+0.2*round(rand(data,p_max)); % Mask
    if abs(sum(preM)) < 1e-5
        maskCheck = 0;
    end
end

% 長めの入力を作成
Jl = zeros(length(preM),dataLen, data); %行がマスク，列が入力uの2次元配列
Jt = zeros(length(preM),dataLen, data);
Jl2 = zeros(length(preM2),dataLen, data); %行がマスク，列が入力uの2次元配列
Jt2 = zeros(length(preM2),dataLen, data);
for step = 1:dataLen
    Jl(:,step) = ul(:,step).*preM;
    Jt(:,step) = ut(:,step).*preM;
    Jl2(:,step) = ul(:,step).*preM2;
    Jt2(:,step) = ut(:,step).*preM2;
end
Jl2=Jl;Jt2=Jt;
% リザーバのデータを毎ステップとって，周期の始点を検出しながらMGの計算
Xl = zeros(1,dataLen*N);
Xt = zeros(1,dataLen*N);
Xl2 = zeros(1,dataLen*N);
Xt2 = zeros(1,dataLen*N);
X_index_l = ones(dataLen,2,data); % 入力した長さを記録
X_index_t = ones(dataLen,2,data);

X_tau = 0; % Xl(:,1)
X_tau2 =0; %Xt(:,1)
X_tau_2 = 0; % Xl(:,1)
X_tau2_2 =0; %Xt(:,1)
u_step_l = 1;
u_step_t = 1;
M_step_l = 1;
M_step_t = 1;
for step = 2:N 
    M_step_l = M_step_l +1;
    dxl = dif(Xl(:,step-1),X_tau, X_tau_2, a, W(1,1), W(1,2), gamma, Jl(M_step_l,u_step_l));
    Xl(:,step) = Xl(:,step-1) + dxl*theta;
    
    dxl2 = dif(Xl2(:,step-1),X_tau_2, X_tau, a, W(2,2), W(2,1), gamma, Jl2(M_step_l,u_step_l));
    Xl2(:,step) = Xl2(:,step-1) + dxl2*theta;
    
    
    
    M_step_t = M_step_t +1;
    dxt = dif(Xt(:,step-1), X_tau2, X_tau2_2, a,  W(1,1), W(1,2), gamma, Jt(M_step_t,u_step_t));
    Xt(:,step) = Xt(:,step-1) + dxt*theta;
    
    dxt2 = dif(Xt2(:,step-1),X_tau2_2, X_tau2, a, W(2,2), W(2,1), gamma, Jt2(M_step_t,u_step_t));
    Xt2(:,step) = Xt2(:,step-1) + dxt2*theta;
end

while u_step_l < dataLen
    
    if step >= p_start(u_step_l+1,1)-1
        X_index_l(u_step_l,2) = M_step_l;
        X_index_l(u_step_l+1, 1) = X_index_l(u_step_l, 1)+M_step_l;
        M_step_l = 0;
        u_step_l =u_step_l+1;
    end
    step = step+1;
    M_step_l = M_step_l +1;
    X_tau = Xl(:,step-N); X_tau_2 = Xl2(:,step-N);
    dxl = dif(Xl(:,step-1),X_tau, X_tau_2, a,  W(1,1), W(1,2), gamma, Jl(M_step_l,u_step_l));
    Xl(:,step) = Xl(:,step-1) + dxl*theta;
    
    dxl2 = dif(Xl2(:,step-1),X_tau_2, X_tau, a, W(2,2), W(2,1), gamma, Jl2(M_step_l,u_step_l));
    Xl2(:,step) = Xl2(:,step-1) + dxl2*theta;
end

step = N;
while u_step_t < dataLen
    if step >= p_start(u_step_t+1,1)-1
        X_index_t(u_step_t,2) = M_step_t;
        X_index_t(u_step_t+1, 1) = X_index_t(u_step_t, 1)+M_step_t;
        M_step_t = 0;
        u_step_t =u_step_t+1;
    end
    step = step+1;
    M_step_t = M_step_t +1;
    X_tau2 = Xt(:,step-N); X_tau2_2 = Xt2(:,step-N);
    dxt = dif(Xt(:,step-1), X_tau2, X_tau2_2, a,  W(1,1), W(1,2), gamma, Jt(M_step_t,u_step_t));
    Xt(:,step) = Xt(:,step-1) + dxt*theta;
    
    dxt2 = dif(Xt2(:,step-1),X_tau2_2, X_tau2, a, W(2,2), W(2,1), gamma, Jt2(M_step_t,u_step_t));
    Xt2(:,step) = Xt2(:,step-1) + dxt2*theta;
end

%% 学習データ
l_interval = N;

x_kl0 = zeros(learnDimension, dataLen);
x_kt0 = zeros(learnDimension, dataLen);
x_kl2 = zeros(learnDimension, dataLen);
x_kt2 = zeros(learnDimension, dataLen);
for k_step = 1:dataLen-1
    x_kl0(:,k_step) = Xl(:,X_index_l(k_step,1):l_interval/learnDimension:X_index_l(k_step,1)+l_interval-1);
    x_kt0(:,k_step) = Xt(:,X_index_t(k_step,1):l_interval/learnDimension:X_index_t(k_step,1)+l_interval-1);
    x_kl2(:,k_step) = Xl2(:,X_index_l(k_step,1):l_interval/learnDimension:X_index_l(k_step,1)+l_interval-1);
    x_kt2(:,k_step) = Xt2(:,X_index_t(k_step,1):l_interval/learnDimension:X_index_t(k_step,1)+l_interval-1);
end

%% ここでvertcat
x_kl = [x_kl0; x_kl2];
x_kt = [x_kt0; x_kt2];

if inputCheck == 0
    x_kl = vertcat(zeros(1,dataLen),x_kl);
    x_kt = vertcat(zeros(1,dataLen),x_kt);
else
    x_kl = vertcat(ul,x_kl);
    x_kt = vertcat(ut,x_kt);
end
if biasCheck == 0
    x_kl = vertcat(zeros(1,dataLen),x_kl);
    x_kt = vertcat(zeros(1,dataLen),x_kt);
else
    x_kl = vertcat(ones(1,dataLen),x_kl);
    x_kt = vertcat(ones(1,dataLen),x_kt);
end
end

%% モデル定義
function dx = dif(X, X_tau, X2_tau,  a, xb, xb2, gamma, J)
dx = -a*X + tanh(xb*X_tau + xb2*X2_tau + gamma*J);
end
