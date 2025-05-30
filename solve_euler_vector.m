% 主程序：从文件读取数据，筛选范围后计算欧拉矢量并输出结果
clear; clc;

% 文件路径
data_file = 'data.txt';  % 输入文件
output_file = 'out.txt';   % 输出文件

% 设置筛选范围（[最小值, 最大值]）
% 以土耳其为例36° 至42° N，经度26° 至45° E
lon_range = [25, 45];  % 经度范围
lat_range = [35, 45];  % 纬度范围

% 统一地球半径 (m)
R = 6371e3;  % 6371 km 转换为 m

% 读取数据文件并处理空格
fid = fopen(data_file, 'r');
if fid == -1
    error('无法打开输入文件!');
end

data = load(data_file);

% 提取数据列
lon = data(:,1);  % 经度
lat = data(:,2);  % 纬度
ve = data(:,3);   % 东向速度分量 (mm/yr)
vn = data(:,4);   % 北向速度分量 (mm/yr)

% 筛选经纬度范围内的数据
in_range = (lon >= lon_range(1) & lon <= lon_range(2) & ...
            lat >= lat_range(1) & lat <= lat_range(2));

lon_filtered = lon(in_range);
lat_filtered = lat(in_range);
ve_filtered = ve(in_range);
vn_filtered = vn(in_range);

% 显示筛选结果
fprintf('原始数据点: %d 个\n', length(lon));
fprintf('筛选后数据点: %d 个\n', length(lon_filtered));
fprintf('筛选范围: 经度 %.2f° 到 %.2f°, 纬度 %.2f° 到 %.2f°\n', ...
    lon_range(1), lon_range(2), lat_range(1), lat_range(2));

% 检查筛选后是否有数据
if length(lon_filtered) < 3
    error('筛选后的数据点数量不足，至少需要3个点来计算欧拉矢量');
end

% 创建stations结构体数组
nstations = length(lon_filtered);
stations = struct('lon', lon_filtered, 'lat', lat_filtered, ...
                  've', ve_filtered, 'vn', vn_filtered);

% 计算欧拉矢量
[omega_true, pole_lat, pole_lon, rate] = euler_vector(stations, R);

% 打印结果
fprintf('欧拉极位置: 经度 = %.4f°, 纬度 = %.4f°\n', pole_lon, pole_lat);
fprintf('旋转速率: %.6f°/Myr\n', rate);

% 将结果写入输出文件
fid = fopen(output_file, 'w');
if fid == -1
    error('无法打开输出文件!');
end

% 写入文件头
fprintf(fid, '# 欧拉矢量计算结果\n');
fprintf(fid, '# 筛选范围: 经度 %.2f° 到 %.2f°, 纬度 %.2f° 到 %.2f°\n', ...
    lon_range(1), lon_range(2), lat_range(1), lat_range(2));
fprintf(fid, '# 欧拉极位置: 经度 = %.4f°, 纬度 = %.4f°\n', pole_lon, pole_lat);
fprintf(fid, '# 旋转速率: %.6f°/Myr\n', rate);
fprintf(fid, '#\n');
fprintf(fid, '# 格式: 经度 纬度 东向速度(mm/yr) 北向速度(mm/yr) 预测东向速度(mm/yr) 预测北向速度(mm/yr) 残差东向(mm/yr) 残差北向(mm/yr)\n');

% 计算每个点的预测速度和残差 (单位转换: m/yr -> mm/yr)
pred_ve = zeros(nstations, 1);
pred_vn = zeros(nstations, 1);
resid_ve = zeros(nstations, 1);
resid_vn = zeros(nstations, 1);

for i = 1:nstations
    lon_rad = deg2rad(stations(i).lon);
    lat_rad = deg2rad(stations(i).lat);
    
    % 使用点乘 .*
    x = cos(lat_rad) .* cos(lon_rad);
    y = cos(lat_rad) .* sin(lon_rad);
    z = sin(lat_rad);
    
    % 计算预测速度 (使用角速度矢量omega_true, 单位: m/yr)
    v_east_m_per_yr = R * (omega_true(3) * y - omega_true(2) * z);
    v_north_m_per_yr = R * (omega_true(1) * z - omega_true(3) * x);
    
    % 转换为mm/yr
    pred_ve(i) = v_east_m_per_yr * 1e3;
    pred_vn(i) = v_north_m_per_yr * 1e3;
    
    % 计算残差 (单位: mm/yr)
    resid_ve(i) = stations(i).ve - pred_ve(i);
    resid_vn(i) = stations(i).vn - pred_vn(i);
end

% 写入数据
for i = 1:nstations
    fprintf(fid, '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        stations(i).lon, stations(i).lat, ...
        stations(i).ve, stations(i).vn, ...
        pred_ve(i), pred_vn(i), ...
        resid_ve(i), resid_vn(i));
end

fclose(fid);
fprintf('结果已写入文件: %s\n', output_file);

% 函数：求解欧拉矢量（普通最小二乘法）
function [omega_true, pole_lat, pole_lon, rate] = euler_vector(stations, R)
    nstations = length(stations);
    b = zeros(2*nstations, 1);
    A = zeros(2*nstations, 3);
    
    for i = 1:nstations
        lon = deg2rad(stations.lon(i));
        lat = deg2rad(stations.lat(i));
        
%         % 使用点乘 .*
%         x = cos(lat) .* cos(lon);
%         y = cos(lat) .* sin(lon);
%         z = sin(lat);
        
        idx = 2*i - 1;
        b(idx) = stations.ve(i);   % 东向速度 (mm/yr)
        b(idx+1) = stations.vn(i); % 北向速度 (mm/yr)
        
        % 设计矩阵
        A(idx,:)   = [R*cos(lon)*sin(lat), R*-sin(lat)*sin(lon), R*cos(lat)];  % 对应: ve = y*Ωz - z*Ωy
        A(idx+1,:) = [R*sin(lon),R*cos(lat) , 0];  % 对应: vn = z*Ωx - x*Ωz
    end
    
    % 求解最小二乘问题
    Omega = (A' * A) \ (A' * b);
    
    % 计算真实的角速度矢量 (rad/year)
    omega_true = Omega / (R * 1e3);  % R单位: m, 转换为mm需要乘以1e3
    
    % 计算旋转速率 (°/Myr)
    omega_mag = norm(omega_true);
    rate = rad2deg(omega_mag) * 1e6;
    
    % 计算欧拉极位置
    if omega_mag > 0
        dir_vec = omega_true / omega_mag;  % 归一化方向向量
    else
        dir_vec = [0, 0, 0];
    end
    
    pole_lat = rad2deg(asin(dir_vec(3)));
    pole_lon = rad2deg(atan2(dir_vec(2), dir_vec(1)));
end