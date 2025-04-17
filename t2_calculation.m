function [] = t2_calculation(S_data)
% FBG三维曲线重建 - 根据FBG传感器数据重建三维曲线形状
% 输入：S_data (包含s, kappa_a, kappa_b, kappa_c列的数据表)

%% 1. 数据准备
% 提取并确保数据为列向量
s_raw = S_data{:,2};          % 弧长s列
kappa_a = S_data{:,3};        % FBG阵列a的曲率
kappa_b = S_data{:,4};        % FBG阵列b的曲率
kappa_c = S_data{:,5};        % FBG阵列c的曲率

% 数据清洗步骤
% (1) 移除NaN值
valid_idx = ~isnan(s_raw) & ~isnan(kappa_a) & ~isnan(kappa_b) & ~isnan(kappa_c);
s_raw = s_raw(valid_idx);
kappa_a = kappa_a(valid_idx);
kappa_b = kappa_b(valid_idx);
kappa_c = kappa_c(valid_idx);

% (2) 确保弧长单调递增
[s_raw, sort_idx] = sort(s_raw);
kappa_a = kappa_a(sort_idx);
kappa_b = kappa_b(sort_idx);
kappa_c = kappa_c(sort_idx);

% (3) 检查并处理重复点
[~, unique_idx] = unique(s_raw);
if length(unique_idx) < length(s_raw)
    warning('发现重复的弧长值，使用平均值合并');
    [s_raw, ~, ic] = unique(s_raw);
    kappa_a = accumarray(ic, kappa_a, [], @mean);
    kappa_b = accumarray(ic, kappa_b, [], @mean);
    kappa_c = accumarray(ic, kappa_c, [], @mean);
end

%% 2. 生成高密度插值点
n_interp = 10; % 插值密度因子
s_fine = linspace(min(s_raw), max(s_raw), n_interp*length(s_raw))';

%% 3. 曲率插值计算
try
    % 计算曲率分量
    kx_raw = kappa_a - 0.5*(kappa_b + kappa_c);
    ky_raw = (sqrt(3)/2)*(kappa_b - kappa_c);
    
    % 使用保形插值(pchip)提高稳定性
    k_x = pchip(s_raw, kx_raw, s_fine);
    k_y = pchip(s_raw, ky_raw, s_fine);
    
    % 处理可能的NaN值
    k_x(isnan(k_x)) = 0;
    k_y(isnan(k_y)) = 0;
    
    % 计算总曲率
    kappa = sqrt(k_x.^2 + k_y.^2);
    kappa(kappa < eps) = eps; % 避免除零
    
catch ME
    error('曲率插值失败: %s', ME.message);
end

%% 4. 法向量计算
N = zeros(length(s_fine), 3);
N(:,1) = k_x ./ kappa;  % X方向分量
N(:,2) = k_y ./ kappa;  % Y方向分量
% Z分量保持为0

%% 5. Frenet标架数值积分（RK4方法）
% 初始化变量
num_points = length(s_fine);
T = zeros(num_points, 3);  % 切向量
B = zeros(num_points, 3);  % 副法向量
r = zeros(num_points, 3);  % 空间坐标

% 初始条件（起点在原点，初始方向沿z轴）
T(1,:) = [0 0 1];          
B(1,:) = cross(T(1,:), N(1,:)); 
r(1,:) = [0 0 0];          

% 计算平均步长
h = mean(diff(s_fine));

% RK4积分主循环
for i = 2:num_points
    %% RK4求解切向量更新
    % 阶段1
    k1_T = h * kappa(i-1) * N(i-1,:);
    
    % 阶段2（中点）
    T_mid = T(i-1,:) + k1_T/2;
    T_mid = T_mid/norm(T_mid);
    k2_T = h * mean(kappa(i-1:i)) * mean(N(i-1:i,:));
    
    % 阶段3（改进中点）
    T_mid = T(i-1,:) + k2_T/2;
    T_mid = T_mid/norm(T_mid);
    k3_T = h * mean(kappa(i-1:i)) * mean(N(i-1:i,:));
    
    % 阶段4（终点）
    T_end = T(i-1,:) + k3_T;
    T_end = T_end/norm(T_end);
    k4_T = h * kappa(i) * N(i,:);
    
    % 加权平均更新切向量
    T(i,:) = T(i-1,:) + (k1_T + 2*k2_T + 2*k3_T + k4_T)/6;
    T(i,:) = T(i,:)/norm(T(i,:));
    
    %% 更新副法向量和位置
    B(i,:) = cross(T(i,:), N(i,:));
    r(i,:) = r(i-1,:) + h * T(i-1,:);
end

%% 6. 结果可视化
figure('Position', [100 100 1200 500]);

% 3D曲线图
subplot(1,2,1);
plot3(r(:,1), r(:,2), r(:,3), 'b-', 'LineWidth', 1.5);
hold on;
scatter3(r(1:100:end,1), r(1:100:end,2), r(1:100:end,3), 30, 'r', 'filled');
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
title('Reconstructed 3D Curve');
axis equal tight; grid on;
view(30,30); rotate3d on;

% 曲率变化图
subplot(1,2,2);
plot(s_fine, kappa, 'LineWidth', 1.5);
xlabel('Arc Length s [mm]'); ylabel('Curvature \kappa [1/mm]');
title('Curvature Profile');
grid on;

%% 7. 结果输出
output_table = array2table([s_fine, r], ...
    'VariableNames', {'ArcLength_mm', 'X_mm', 'Y_mm', 'Z_mm'});
writetable(output_table, 'Reconstructed_Curve.csv');
disp('First 5 reconstructed points:');
disp(head(output_table, 5));
end