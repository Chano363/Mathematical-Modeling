%% 1. 数据预处理
s = data{:,1};               % 提取弧长列
kappa_data = S_data{:,2:n+1};  % 提取n列曲率数据

% 数据有效性检查
assert(size(kappa_data,2) == n, '曲率数据列数与指定的n值不符');
assert(n >= 3, 'FBG阵列数量n必须≥3');

% 数据清洗
valid_idx = all(~isnan([s, kappa_data]), 2);
s = s(valid_idx);
kappa_data = kappa_data(valid_idx,:);

% 按弧长排序
[s, sort_idx] = sort(s);
kappa_data = kappa_data(sort_idx,:);

% 处理重复点
[~, unique_idx] = unique(s);
if length(unique_idx) < length(s)
    [s, ~, ic] = unique(s);
    kappa_data = accumarray(ic, (1:length(ic))', [], @(x) mean(kappa_data(x,:),1));
end

%% 2. 高密度插值
n_interp = 10; % 插值密度因子
s_fine = linspace(min(s), max(s), n_interp*length(s))';

% 对各FBG通道分别插值
kappa_interp = zeros(length(s_fine), n);
for i = 1:n
    kappa_interp(:,i) = pchip(s, kappa_data(:,i), s_fine);
end
kappa_interp(isnan(kappa_interp)) = 0;

%% 3. 方向矩阵构建
phi = linspace(0, 360, n+1)'; phi(end) = []; % 等间隔角度分布
D = [cosd(phi), sind(phi), zeros(n,1)];       % n×3方向矩阵

%% 4. 曲率向量估计（最小二乘法）
kappa_vectors = zeros(length(s_fine), 3);
for i = 1:length(s_fine)
    k = kappa_interp(i,:)';
    kappa_vectors(i,:) = (D'*D) \ (D'*k); % 最小二乘解
end

% 计算曲率大小
kappa = sqrt(sum(kappa_vectors.^2, 2));
kappa(kappa < eps) = eps; % 避免除零

%% 5. 单位法向量计算
N = kappa_vectors ./ kappa;
N(isnan(N)) = 0; % 处理零曲率情况

%% 6. RK4积分重建曲线
num_points = length(s_fine);
T = zeros(num_points, 3);  % 切向量
B = zeros(num_points, 3);  % 副法向量
r = zeros(num_points, 3);  % 空间坐标

% 初始条件（起点在原点，初始方向沿z轴）
T(1,:) = [0 0 1];
B(1,:) = cross(T(1,:), N(1,:));
r(1,:) = [0 0 0];

h = mean(diff(s_fine)); % 平均步长

% RK4积分主循环
for i = 2:num_points
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
    
    % 加权平均更新
    T(i,:) = T(i-1,:) + (k1_T + 2*k2_T + 2*k3_T + k4_T)/6;
    T(i,:) = T(i,:)/norm(T(i,:));
    
    % 更新副法向量和位置
    B(i,:) = cross(T(i,:), N(i,:));
    r(i,:) = r(i-1,:) + h * T(i-1,:);
end

%% 7. 结果可视化
figure('Position', [100 100 1200 500]);

% 3D曲线
subplot(1,2,1);
plot3(r(:,1), r(:,2), r(:,3), 'b-', 'LineWidth', 1.5);
hold on;
scatter3(r(1:100:end,1), r(1:100:end,2), r(1:100:end,3), 30, 'r', 'filled');
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
title(['Reconstructed 3D Curve (n=', num2str(n), ' FBG arrays)']);
axis equal tight; grid on;
view(30,30); rotate3d on;

% 曲率变化
subplot(1,2,2);
plot(s_fine, kappa, 'LineWidth', 1.5);
xlabel('Arc Length [mm]'); ylabel('Curvature [1/mm]');
title('Curvature Profile');
grid on;

%% 8. 结果输出
output_table = array2table([s_fine, r], ...
    'VariableNames', {'ArcLength_mm', 'X_mm', 'Y_mm', 'Z_mm'});
writetable(output_table, ['Reconstructed_Curve_n', num2str(n), '.csv']);
disp(['First 5 reconstructed points (n=', num2str(n), ' arrays):']);
disp(head(output_table, 5));

%% 9. 条件数分析（可选）
cond_number = cond(D'*D);
fprintf('\n方向矩阵条件数: %.2f\n', cond_number);
if cond_number > 1e4
    warning('高条件数可能导致数值不稳定');
end
end