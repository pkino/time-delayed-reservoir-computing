function  [step_half, N, l_num, x_kl, x_kt, ul, ut, l_start] = Res(a, b, c1, gamma)
%% 初期設定
% data*step_allの行列のデータを入力
data = 1; % 1ステップあたりの入力データ数
step_half = 1500; % 初期値は3200 論文執筆時は実は1500
step_all = 2*step_half; % データ長=ステップ数
theta = 0.01; % [s] 0.2 0.01

tau = 80; %[s] 400*theta
N = tau/theta; % delayのstep数
l_num = 400;
% lin = 1;
% eta = 0.5; % 1.8 2 1.3
% c1 =-0.1; %0=線形モデル、-0.1=非線形モデル
% gamma =0.5;

%データの長さ
step_half =step_half*N;
step_all = step_all*N;




%% 目標データの生成
data_length = step_half/N+200; % 少し長めに準備

% 格納変数・初期値
ul = 0.5*rand(data,data_length); % 入力データを学習とテストで作り変える
ut = 0.5*rand(data,data_length);

% MG = load('MackeyGlass_t17.txt');
% ul = MG(1:data_length)';
% ut = MG(data_length+1:data_length*2)';


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
    dxl = dif(Xl(:,step-1),X_tau, a, b, c1, gamma, Jl(M_step_l,u_step_l));
    Xl(:,step) = Xl(:,step-1) + dxl*theta;
    
    M_step_t = M_step_t +1;
    dxt = dif(Xt(:,step-1),X_tau2, a, b, c1, gamma, Jt(M_step_t,u_step_t));
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
    dxl = dif(Xl(:,step-1), X_tau, a, b, c1, gamma, Jl(M_step_l,u_step_l));
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
    dxt = dif(Xt(:,step-1), X_tau2, a, b, c1, gamma, Jt(M_step_t,u_step_t));
    Xt(:,step) = Xt(:,step-1) + dxt*theta;
end

%% 学習データから重みを求める memoryFunction
l_interval = N;
l_start = 100 + 50; % X(X0)が安定したところから学習スタート
x_kl = zeros(l_num, data_length);
x_kt = zeros(l_num, data_length);
for k_step = 1:data_length-1
    x_kl(:,k_step) = Xl(:,X_index_l(k_step,1):l_interval/l_num:X_index_l(k_step,1)+l_interval-1);
    x_kt(:,k_step) = Xt(:,X_index_t(k_step,1):l_interval/l_num:X_index_t(k_step,1)+l_interval-1);
end

%% テスト
% テストまで行う場合はこのコメントアウトを外す
% taskDelay = 50;
% RC_MF_pinv(taskDelay, step_half, N, l_num,x_kl, x_kt, ul, ut, l_start)

end

%% モデル定義
function dx = dif(X, X_tau, a, b, c1, gamma, J)
dx = -a*X + b*X_tau + c1*X^3 + gamma*J;
end
