function reconstruct_curve(k_data, s, angles)
    % 提取光纤间夹角
    alpha = angles(1);
    beta = angles(2);
    % 获取检测点数量
    num_points = length(s);
    % 初始化曲率、角度和挠率数组
    curvature = zeros(num_points, 1);
    angle = zeros(num_points, 1);
    torsion = zeros(num_points, 1);

    % 高精度曲率转换
    for i = 1:num_points
        % 调用函数计算二维等效曲率分量
        [curvature_x, curvature_y] = corrected_curvature(k_data(i, 1), k_data(i, 2), k_data(i, 3), alpha, beta);
        % 计算合成后的曲率大小
        curvature(i) = sqrt(curvature_x^2 + curvature_y^2);
        % 计算曲率方向角度
        angle(i) = atan2(curvature_y, curvature_x);
    end

    % 中心差分法计算挠率
    for i = 2:num_points - 1
        % 根据相邻点的角度差和弧长差计算挠率
        torsion(i) = (angle(i + 1) - angle(i - 1)) / (s(i + 1) - s(i - 1));
    end
    % 计算第一个点的挠率（向前差分）
    torsion(1) = (angle(2) - angle(1)) / (s(2) - s(1));
    % 计算最后一个点的挠率（向后差分）
    torsion(num_points) = (angle(num_points) - angle(num_points - 1)) / (s(num_points) - s(num_points - 1));

    % 初始化切向量、法向量和副法向量
    tangent_vector = [0; 0; 1];
    normal_vector = [cos(angle(1)); sin(angle(1)); 0];
    binormal_vector = cross(tangent_vector, normal_vector);
    % 初始化曲线坐标矩阵
    r = zeros(num_points, 3);

    % 使用RK4积分Frenet-Serret公式计算曲线坐标
    for i = 1:num_points - 1
        % 计算相邻检测点间的弧长增量
        ds = s(i + 1) - s(i);
        % 根据切向量和弧长增量更新曲线坐标
        r(i + 1, :) = r(i, :) + (tangent_vector * ds)';
        % 更新切向量、法向量和副法向量
        [tangent_vector, normal_vector, binormal_vector] = update_frenet_frame(tangent_vector, normal_vector, binormal_vector, curvature(i), torsion(i), ds);
    end

    % 绘制空间曲线
    figure;
    hold on;
    scatter3(r(:,1), r(:,2), r(:,3), 'ro', 'DisplayName', 'Reconstructed Curve Points');
    xlabel('X - Coordinate');
    ylabel('Y - Coordinate');
    zlabel('Z - Coordinate');
    title('Reconstructed Space Curve');
    legend;
    grid on;
    hold off;
end