function [] = t3_calculation(S_data)
%% 1. 数据预处理
s = S_data{:,2};          % 提取弧长数据列
kappa_a = S_data{:,3};    % 提取FBG阵列a的曲率
kappa_b = S_data{:,4};    % 提取FBG阵列b的曲率
kappa_c = S_data{:,5};    % 提取FBG阵列c的曲率

% 数据清洗：去除无效数据点
valid_idx = ~isnan(s) & ~isnan(kappa_a) & ~isnan(kappa_b) & ~isnan(kappa_c);
[s, sort_idx] = sort(s(valid_idx));
kappa_a = kappa_a(valid_idx);
kappa_b = kappa_b(valid_idx);
kappa_c = kappa_c(valid_idx);
kappa_a = kappa_a(sort_idx);
kappa_b = kappa_b(sort_idx);
kappa_c = kappa_c(sort_idx);

% 处理重复数据点：取平均值
[~, unique_idx] = unique(s);
if length(unique_idx) < length(s)
    [s, ~, ic] = unique(s);
    kappa_a = accumarray(ic, kappa_a, [], @mean);
    kappa_b = accumarray(ic, kappa_b, [], @mean);
    kappa_c = accumarray(ic, kappa_c, [], @mean);
end

%% 2. 高密度插值处理
n_interp = 10; % 插值密度系数
s_fine = linspace(min(s), max(s), n_interp*length(s))'; % 生成精细弧长网格

%% 3. 曲率计算与处理
% 计算曲率分量
kx_raw = kappa_a - 0.5*(kappa_b + kappa_c); % x方向曲率分量
ky_raw = (sqrt(3)/2)*(kappa_b - kappa_c);   % y方向曲率分量

% 使用保形插值(pchip)提高稳定性
k_x = pchip(s, kx_raw, s_fine); % x方向曲率插值
k_y = pchip(s, ky_raw, s_fine); % y方向曲率插值

% 处理可能的NaN值
k_x(isnan(k_x)) = 0;
k_y(isnan(k_y)) = 0;

% 计算总曲率大小
kappa = sqrt(k_x.^2 + k_y.^2);
kappa(kappa < eps) = eps; % 避免除零错误

%% 4. 坐标系计算
% 计算单位法向量
N = zeros(length(s_fine), 3);
N(:,1) = k_x ./ kappa; % x分量
N(:,2) = k_y ./ kappa; % y分量
% z分量保持为0

%% 5. 三维曲线重建（RK4积分法）
num_points = length(s_fine);
T = zeros(num_points, 3); % 切向量
B = zeros(num_points, 3); % 副法向量
r = zeros(num_points, 3); % 空间坐标

% 初始条件：起点在原点，初始方向沿z轴
T(1,:) = [0 0 1];          
B(1,:) = cross(T(1,:), N(1,:)); 
r(1,:) = [0 0 0];          

% 计算平均步长
h = mean(diff(s_fine));

% RK4积分主循环
for i = 2:num_points
    %% RK4法更新切向量
    % 第一阶段
    k1_T = h * kappa(i-1) * N(i-1,:);
    
    % 第二阶段（中点）
    T_mid = T(i-1,:) + k1_T/2;
    T_mid = T_mid/norm(T_mid);
    k2_T = h * mean(kappa(i-1:i)) * mean(N(i-1:i,:));
    
    % 第三阶段（改进中点）
    T_mid = T(i-1,:) + k2_T/2;
    T_mid = T_mid/norm(T_mid);
    k3_T = h * mean(kappa(i-1:i)) * mean(N(i-1:i,:));
    
    % 第四阶段（终点）
    T_end = T(i-1,:) + k3_T;
    T_end = T_end/norm(T_end);
    k4_T = h * kappa(i) * N(i,:);
    
    % 加权平均更新切向量
    T(i,:) = T(i-1,:) + (k1_T + 2*k2_T + 2*k3_T + k4_T)/6;
    T(i,:) = T(i,:)/norm(T(i,:));
    
    %% 更新副法向量和位置坐标
    B(i,:) = cross(T(i,:), N(i,:));
    r(i,:) = r(i-1,:) + h * T(i-1,:);
end

%% 6. 结果可视化
figure('Position', [100 100 1200 500]);

% 三维曲线图
subplot(1,2,1);
plot3(r(:,1), r(:,2), r(:,3), 'b-', 'LineWidth', 1.5);
hold on;
scatter3(r(1:100:end,1), r(1:100:end,2), r(1:100:end,3), 30, 'r', 'filled');
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
title('Reconstructed 3D Shape');
axis equal tight; grid on;
view(30,30); rotate3d on;

% 曲率分布图
subplot(1,2,2);
plot(s_fine, kappa, 'LineWidth', 1.5);
xlabel('Arc Length [mm]'); ylabel('Curvature [1/mm]');
title('Curvature Profile');
grid on;

%% 7. 结果输出
output_table = array2table([s_fine, r], ...
    'VariableNames', {'ArcLength_mm', 'X_mm', 'Y_mm', 'Z_mm'});
writetable(output_table, 'Reconstructed_Curve_3D.csv');
disp('First 5 reconstructed points:');
disp(head(output_table, 5));
end