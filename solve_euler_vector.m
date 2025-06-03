% 主程序：从文件读取数据，筛选范围后计算欧拉矢量并输出结果
clear; clc;

% 文件路径
data_file = 'data.txt';  % 输入文件
output_file = 'out1.txt'; % 输出文件

% 设置筛选范围（[最小值, 最大值]）
lon_range = [20, 30];  % 经度范围20-50
lat_range = [35, 45];  % 纬度范围35-45

% 统一地球半径 (m)
R = 6370e6;  % 6371 km 转换为 mm

% 读取数据
data = readData(data_file);

% 提取数据列
lon = data(:,1);  % 经度
lat = data(:,2);  % 纬度
ve = data(:,3);   % 东向速度分量 (mm/yr)
vn = data(:,4);   % 北向速度分量 (mm/yr)

% 筛选经纬度范围内的数据
[lon_filtered, lat_filtered, ve_filtered, vn_filtered] = filterData(lon, lat, ve, vn, lon_range, lat_range);

% 显示筛选结果
fprintf('原始数据点: %d 个\n', length(lon));
fprintf('筛选后数据点: %d 个\n', length(lon_filtered));
fprintf('筛选范围: 经度 %.2f° 到 %.2f°, 纬度 %.2f° 到 %.2f°\n', ...
    lon_range(1), lon_range(2), lat_range(1), lat_range(2));

% 检查筛选后是否有足够数据
validateDataLength(lon_filtered);

% 创建stations结构体数组
stations = createStationStruct(lon_filtered, lat_filtered, ve_filtered, vn_filtered);

% 计算欧拉矢量
[omega_true, pole_lat, pole_lon, rate] = euler_vector(stations, R);

% 计算预测速度和残差
[pred_ve, pred_vn, resid_ve, resid_vn] = calculatePredictions(stations, omega_true, R);

% 计算残差统计信息
[resid_stats] = calculateResidualStats(resid_ve, resid_vn);

% 打印结果
printResults(pole_lon, pole_lat, rate, resid_stats);

% 将结果写入输出文件
writeResults(output_file, lon_range, lat_range, pole_lon, pole_lat, rate, stations, pred_ve, pred_vn, resid_ve, resid_vn);

% 可视化结果（移除原始数据展示）
visualizeResults(lon_filtered, lat_filtered, ve_filtered, vn_filtered, lon_range, lat_range, pred_ve, pred_vn, resid_ve, resid_vn);

% 函数：读取数据文件
function data = readData(file_path)
    % 检查文件是否存在
    if ~exist(file_path, 'file')
        error('无法打开输入文件: %s', file_path);
    end
    
    % 读取数据
    try
        data = load(file_path);
    catch e
        error('读取文件失败: %s', e.message);
    end
    
    % 验证数据格式
    if size(data, 2) < 4
        error('数据文件格式错误: 至少需要4列数据 (经度, 纬度, 东向速度, 北向速度)');
    end
end

% 函数：筛选数据
function [lon_filtered, lat_filtered, ve_filtered, vn_filtered] = filterData(lon, lat, ve, vn, lon_range, lat_range)
    in_range = (lon >= lon_range(1) & lon <= lon_range(2) & ...
                lat >= lat_range(1) & lat <= lat_range(2));
    
    lon_filtered = lon(in_range);
    lat_filtered = lat(in_range);
    ve_filtered = ve(in_range);
    vn_filtered = vn(in_range);
end

% 函数：验证数据长度
function validateDataLength(data)
    if length(data) < 3
        error('筛选后的数据点数量不足，至少需要3个点来计算欧拉矢量');
    end
end

% 函数：创建站点结构体
function stations = createStationStruct(lon, lat, ve, vn)
    nstations = length(lon);
    stations = struct('lon', lon, 'lat', lat, 've', ve, 'vn', vn);
end

% 函数：计算预测速度和残差
function [pred_ve, pred_vn, resid_ve, resid_vn] = calculatePredictions(stations, omega_true, R)
    nstations = length(stations.lon);
    A = zeros(2*nstations, 3);
    pred_ve = zeros(nstations,1);
    pred_vn = zeros(nstations, 1);
    resid_ve = zeros(nstations, 1);
    resid_vn = zeros(nstations, 1);
    
    for i = 1:nstations
        lon = deg2rad(stations.lon(i));
        lat = deg2rad(stations.lat(i));
        
        idx = 2*i-1;
        % 设计矩阵保持原始形式
        A(idx,:)   = [-R*cos(lon)*sin(lat), -R*sin(lon)*sin(lat), R*cos(lat)]; 
        A(idx+1,:) = [R*sin(lon), -R*cos(lon), 0];  
        
        % 计算预测速度 (使用角速度矢量omega_true, 单位: rad/yr)
        v = A(2*i-1:2*i,:) * omega_true;
        
        % 转换为mm/yr
        pred_ve(i) = v(1) ;  % 从m/yr转换为mm/yr
        pred_vn(i) = v(2) ;  % 从m/yr转换为mm/yr
        
        % 计算残差 (单位: mm/yr)
        resid_ve(i) = stations.ve(i) - pred_ve(i);
        resid_vn(i) = stations.vn(i) - pred_vn(i);
    end
end

% 函数：计算残差统计信息
function stats = calculateResidualStats(resid_ve, resid_vn)
    % 计算总残差
    resid_total = sqrt(resid_ve.^2 + resid_vn.^2);
    
    % 计算统计量
    stats.mean_ve = mean(resid_ve);
    stats.mean_vn = mean(resid_vn);
    stats.std_ve = std(resid_ve);
    stats.std_vn = std(resid_vn);
    stats.rms = sqrt(mean(resid_total.^2));
    stats.max_resid = max(resid_total);
end

% 函数：打印结果
function printResults(pole_lon, pole_lat, rate, resid_stats)
    fprintf('欧拉极位置: 经度 = %.4f°, 纬度 = %.4f°\n', pole_lon, pole_lat);
    fprintf('旋转速率: %.6f°/Myr\n', rate);
    fprintf('残差统计:\n');
    fprintf('  东向速度: 均值 = %.4f mm/yr, 标准差 = %.4f mm/yr\n', resid_stats.mean_ve, resid_stats.std_ve);
    fprintf('  北向速度: 均值 = %.4f mm/yr, 标准差 = %.4f mm/yr\n', resid_stats.mean_vn, resid_stats.std_vn);
    fprintf('  总残差: RMS = %.4f mm/yr, 最大值 = %.4f mm/yr\n', resid_stats.rms, resid_stats.max_resid);
end

% 函数：写入结果到文件
function writeResults(output_file, lon_range, lat_range, pole_lon, pole_lat, rate, stations, pred_ve, pred_vn, resid_ve, resid_vn)
    fid = fopen(output_file, 'w');
    if fid == -1
        error('无法打开输出文件: %s', output_file);
    end
    
    % 写入文件头
    fprintf(fid, '# 欧拉矢量计算结果\n');
    fprintf(fid, '# 筛选范围: 经度 %.2f° 到 %.2f°, 纬度 %.2f° 到 %.2f°\n', ...
        lon_range(1), lon_range(2), lat_range(1), lat_range(2));
    fprintf(fid, '# 欧拉极位置: 经度 = %.4f°, 纬度 = %.4f°\n', pole_lon, pole_lat);
    fprintf(fid, '# 旋转速率: %.6f°/Myr\n', rate);
    fprintf(fid, '#\n');
    fprintf(fid, '# 格式: 经度 纬度 东向速度(mm/yr) 北向速度(mm/yr) 预测东向速度(mm/yr) 预测北向速度(mm/yr) 残差东向(mm/yr) 残差北向(mm/yr)\n');
    
    % 写入数据
    for i = 1:length(stations.lon)
        fprintf(fid, '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
            stations.lon(i), stations.lat(i), ...
            stations.ve(i), stations.vn(i), ...
            pred_ve(i), pred_vn(i), ...
            resid_ve(i), resid_vn(i));
    end
    
    fclose(fid);
    fprintf('结果已写入文件: %s\n', output_file);
end

% 函数：可视化结果（修改后仅显示筛选后数据）
function visualizeResults(lon_filtered, lat_filtered, ve_filtered, vn_filtered, lon_range, lat_range, pred_ve, pred_vn, resid_ve, resid_vn)
    figure('Position', [100, 100, 1200, 800]);
    
    % 绘制筛选后的数据和速度矢量
    subplot(2, 2, 1);
    quiver(lon_filtered, lat_filtered, ve_filtered, vn_filtered, 0, 'r', 'LineWidth', 1.5);
    title('筛选后的数据和速度矢量');
    xlabel('经度 (°)');
    ylabel('纬度 (°)');
    grid on;
    axis equal;
    
    % 绘制观测速度与预测速度对比
    subplot(2, 2, 2);
    quiver(lon_filtered, lat_filtered, ve_filtered, vn_filtered, 0, 'r', 'LineWidth', 1.5);
    hold on;
    quiver(lon_filtered, lat_filtered, pred_ve, pred_vn, 0, 'g', 'LineWidth', 1.5);
    title('观测速度(红)与预测速度(绿)对比');
    xlabel('经度 (°)');
    ylabel('纬度 (°)');
    grid on;
    axis equal;
    legend('观测速度', '预测速度');
    
    % 绘制残差分布
    subplot(2, 2, 3);
    quiver(lon_filtered, lat_filtered, resid_ve, resid_vn, 0, 'b', 'LineWidth', 1.5);
    title('残差分布');
    xlabel('经度 (°)');
    ylabel('纬度 (°)');
    grid on;
    axis equal;
    
    % 绘制残差统计直方图
    subplot(2, 2, 4);
    resid_total = sqrt(resid_ve.^2 + resid_vn.^2);
    histogram(resid_total, 10, 'Normalization', 'pdf');
    title('残差分布直方图');
    xlabel('残差大小 (mm/yr)');
    ylabel('概率密度');
    grid on;
    
    % 调整子图布局
    sgtitle('板块运动分析结果可视化');
    tight_layout;
end

% 函数：求解欧拉矢量（普通最小二乘法，保持设计矩阵不变）
function [omega_true, pole_lat, pole_lon, rate] = euler_vector(stations, R)
    nstations = length(stations.lon);
    b = zeros(2*nstations, 1);
    A = zeros(2*nstations, 3);
    
    for i = 1:nstations
        lon = deg2rad(stations.lon(i));
        lat = deg2rad(stations.lat(i));
        
        idx = 2*i - 1;
        b(idx) = stations.ve(i);   % 东向速度 (mm/yr)
        b(idx+1) = stations.vn(i); % 北向速度 (mm/yr)
        
        % 保持原始设计矩阵不变
        A(idx,:)   = [-R*cos(lon)*sin(lat), -R*sin(lon)*sin(lat), R*cos(lat)]; 
        A(idx+1,:) = [R*sin(lon), -R*cos(lon), 0];  
    end
    
    % 求解最小二乘问题
    Omega = (A' * A) \ (A' * b);
    
    % 计算真实的角速度矢量 (rad/year)
    omega_true = Omega;
    
    % 计算欧拉极位置（从角速度矢量转换）
    omega_mag = norm(omega_true);
    
    if omega_mag == 0
        pole_lat = 0;
        pole_lon = 0;
    else
        % 从角速度矢量计算欧拉极
        % 注意：这里需要根据设计矩阵的具体形式调整计算方式
        % 此处采用标准转换方法
        pole_lat_rad = asin(omega_true(3) / omega_mag);
        pole_lon_rad = atan2(omega_true(2), omega_true(1));
        
        % 转换为度
        pole_lat = rad2deg(pole_lat_rad);
        pole_lon = rad2deg(pole_lon_rad);
        
        % 调整经度范围到 [-180, 180]
        if pole_lon > 180
            pole_lon = pole_lon - 360;
        elseif pole_lon < -180
            pole_lon = pole_lon + 360;
        end
    end
    
    % 计算旋转速率 (°/Myr)
    rate = omega_mag * (180 / pi) * 1e6;
end
