function t4_calculation(data, n)
%% 输入参数：
% data: 表格数据，第1列为弧长s，后续列为各FBG曲率
% n: FBG阵列数量（n≥3）

%% 1. 数据预处理
s = data{:,1};          % 提取弧长
kappa_data = data{:,2:end}; % 提取曲率数据

% 验证数据维度
assert(size(kappa_data,2) == n, 'Curvature data columns mismatch specified n');

% 数据清洗
valid_idx = all(~isnan([s, kappa_data]),2);
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

%% 2. 方向矩阵生成（等角度分布）
theta = linspace(0, 360, n+1)'; theta(end) = []; % 等分圆周
D = [cosd(theta), sind(theta), zeros(n,1)];      % 生成n×3方向矩阵

% 可视化方向分布
figure
polarplot(deg2rad(theta), ones(n,1), 'ro-', 'LineWidth',2)
title(['FBG Array Orientation (n=',num2str(n),')'])

%% 3. 曲率向量估计（正则化最小二乘）
lambda = 1e-6; % 正则化系数
kappa_vectors = zeros(length(s),3);
for i = 1:length(s)
    k = kappa_data(i,:)';
    kappa_vectors(i,:) = (D'*D + lambda*eye(3)) \ (D'*k);
end

%% 4. 三维曲线重建（RK4积分）
T = zeros(length(s),3); % 切向量
N = kappa_vectors ./ vecnorm(kappa_vectors,2,2); % 法向量
r = zeros(length(s),3); % 空间坐标

% 初始条件
T(1,:) = [0 0 1];      % 初始沿Z轴方向
r(1,:) = [0 0 0];      % 起点在原点

h = mean(diff(s)); % 平均步长

% 主积分循环
for i = 2:length(s)
    % RK4积分参数
    k1 = h * vecnorm(kappa_vectors(i-1,:)) * N(i-1,:);
    
    T_mid = T(i-1,:) + k1/2;
    T_mid = T_mid / norm(T_mid);
    k2 = h * mean(vecnorm(kappa_vectors(i-1:i,:),2,2)) * mean(N(i-1:i,:));
    
    T_end = T(i-1,:) + k2;
    T_end = T_end / norm(T_end);
    k3 = h * vecnorm(kappa_vectors(i,:)) * N(i,:);
    
    % 更新切向量
    T(i,:) = T(i-1,:) + (k1 + 2*k2 + k3)/6;
    T(i,:) = T(i,:) / norm(T(i,:));
    
    % 更新位置
    r(i,:) = r(i-1,:) + h * T(i-1,:);
end

%% 5. 结果可视化
figure('Position', [100 100 1200 500])

% 三维曲线
subplot(1,2,1)
plot3(r(:,1), r(:,2), r(:,3), 'b-', 'LineWidth',1.5)
hold on
quiver3(r(1:50:end,1), r(1:50:end,2), r(1:50:end,3),...
        T(1:50:end,1), T(1:50:end,2), T(1:50:end,3),...
        0.5, 'Color',[0.2 0.6 0.2])
xlabel('X [mm]'), ylabel('Y [mm]'), zlabel('Z [mm]')
title(['3D Reconstructed Curve (n=',num2str(n),' FBG)'])
axis equal tight
grid on
view(-30,30)

% 曲率分量
subplot(1,2,2)
hold on
plot(s, kappa_vectors(:,1), 'r-', 'LineWidth',1.5)
plot(s, kappa_vectors(:,2), 'g-', 'LineWidth',1.5)
plot(s, kappa_vectors(:,3), 'b-', 'LineWidth',1.5)
xlabel('Arc Length [mm]'), ylabel('Curvature [1/mm]')
title('Triaxial Curvature Components')
legend({'X','Y','Z'}, 'Location','best')
grid on

%% 6. 重建质量评估
residuals = zeros(length(s),1);
for i = 1:length(s)
    k_pred = D * kappa_vectors(i,:)';
    residuals(i) = norm(kappa_data(i,:)' - k_pred);
end

fprintf('\n===== Reconstruction Report =====\n')
fprintf('Mean Residual: %.4f\n', mean(residuals))
fprintf('Max Residual: %.4f\n', max(residuals))
fprintf('Position STD: X=%.3f, Y=%.3f, Z=%.3f\n',...
    std(r(:,1)), std(r(:,2)), std(r(:,3)))

%% 7. 数据输出
output = table(s, r(:,1), r(:,2), r(:,3),...
    'VariableNames', {'ArcLength','X','Y','Z'});
writetable(output, ['Reconstructed_n',num2str(n),'.csv'])
end