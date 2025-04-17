clc;
clear;
%% 加载并预处理数据
data_from_excel = readtable('data.xlsx','sheet',"问题1");
s_data = data_from_excel{:, 2};    % 提取弧长s列数据
kappa_data = data_from_excel{:, 3}; % 提取曲率kappa列数据

% 数据有效性检查
assert(all(diff(s_data) > 0), '弧长s必须严格递增！');
assert(~any(isnan(kappa_data)), '曲率数据包含NaN！');

%% 曲率归一化处理（避免数值过大）
kappa_data = kappa_data / max(abs(kappa_data)); % 将曲率归一化到[-1, 1]范围

%% 创建曲率插值函数（使用保形插值'pchip'方法）
kappa_interp = @(s) interp1(s_data, kappa_data, s, 'pchip', 'extrap');

%% 定义微分方程组
ode_func = @(s, y) [
    kappa_interp(s);   % 第一个方程：dθ/ds = κ(s)
    cos(y(1));         % 第二个方程：dx/ds = cosθ 
    sin(y(1))          % 第三个方程：dy/ds = sinθ
];

% 设置初始条件：[初始角度θ0, 初始x坐标, 初始y坐标]
initial_conditions = [0; 0; 0]; % θ=0表示初始方向水平向右

%% 求解微分方程（使用ode45自动步长算法）
[s_sol, y_sol] = ode45(ode_func, [0, max(s_data)], initial_conditions);

% 提取计算结果
theta_sol = y_sol(:, 1); % 角度θ的解
x_sol = y_sol(:, 2);     % x坐标的解
y_sol = y_sol(:, 3);     % y坐标的解

%% 绘制重建的曲线形状（英文输出）
figure;
plot(x_sol, y_sol, 'LineWidth', 1.5, 'Color', [0, 0.5, 0.8]);
hold on;
scatter(x_sol(1:50:end), y_sol(1:50:end), 30, 'r', 'filled');
axis equal; % 保持x和y轴比例一致
xlabel('x axis [mm]');
ylabel('y axis [mm]');
title('The curve shape reconstructed based on curvature');
grid on;
legend('Reconstructed curve', 'Sampling points', 'Location', 'best');

%% （可选）绘制角度θ随弧长的变化曲线
% figure;
% plot(s_sol, theta_sol, 'LineWidth', 1.5);
% xlabel('Arc length s [mm]');
% ylabel('Angle θ [rad]');
% title('Angle variation along arc length');
% grid on;