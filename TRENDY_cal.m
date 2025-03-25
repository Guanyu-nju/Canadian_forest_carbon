clc
clear
%% TRENDY GPP

% Define list of file names (or dynamically fetch file list if applicable)
file_names = {
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CABLE-POP_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLASSIC_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLM6.0_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\EDv3_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\E3SM_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\IBIS_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ISBA-CTRIP_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\JSBACH_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPJmL_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPX-Bern_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\OCN_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ORCHIDEE_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\SDGVM_S3_gpp.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\VISIT_S3_gpp.nc",  
    };

% 定义时间范围
years = 2015:2023; % 提取的年份范围
num_years = length(years);
num_months = 12;

% 动态读取第一个文件的维度
gpp_sample = ncread(file_names{1}, 'gpp');
[lon_size, lat_size, ~] = size(gpp_sample);

% 初始化矩阵，存储月数据和年数据
annual_GPP_combined = nan(180, 360, num_months, num_years, length(file_names)); % 月数据
annual_GPP_sum = nan(180, 360, num_years, length(file_names));                 % 年数据

% 定义每个月的天数（平年为基准）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% 遍历每个文件
for f = 1:length(file_names)
    f
    % 读取变量 'gpp'
    gpp_data = ncread(file_names{f}, 'gpp'); % 假设 gpp 数据是 3D 的：lon x lat x time
    [~, ~, time_dim] = size(gpp_data);

    % 特殊处理 ISBA 文件
    if contains(file_names{f}, 'ISBA')
        % ISBA 文件：360x150 -> 150x360，添加 NaN 补齐为 180x360
        gpp_data_rotated = rot90(gpp_data, 1); % 旋转为 150x360
        gpp_data_rotated = cat(1, gpp_data_rotated, nan(30, 360, time_dim)); % 补齐为 180x360
    elseif contains(file_names{f}, 'CLM6')
        gpp_data_rotated = rot90(gpp_data, 1);
        gpp_data_rotated=[gpp_data_rotated(:,181:360,:),gpp_data_rotated(:,1:180,:)];
    else
        % 普通文件：360x180 -> 180x360
        gpp_data_rotated = rot90(gpp_data, 1);
    end

    % 确定提取的时间范围
    start_idx = time_dim - num_years * num_months + 1; % 从 2015 年开始
    end_idx = time_dim;                               % 到 2023 年最后一个月

    % 提取 2015-2023 年的数据
    gpp_sub_data = gpp_data_rotated(:, :, start_idx:end_idx); % 提取最后 9 年的数据

    % 处理数据（按月累计并考虑闰年）
    for y = 1:num_years
        current_year = years(y);

        % 检查是否为闰年
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % 闰年 2 月有 29 天
        else
            days_in_month(2) = 28; % 平年 2 月有 28 天
        end

        % 计算当年的每月秒数
        seconds_per_month = days_in_month * 24 * 3600;

        % 提取当年数据
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        gpp_yearly_data = gpp_sub_data(:, :, month_start_idx:month_end_idx);

        % 转换单位并累计
        gpp_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            gpp_monthly(:, :, m) = gpp_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % 存入月度矩阵
        annual_GPP_combined(:, :, :, y, f) = gpp_monthly;

        % 将 12 个月的数据求和，存入年度矩阵
        annual_GPP_sum(:, :, y, f) = sum(gpp_monthly, 3); % 对第 3 维求和
    end
end


%% TRENDY Ra

% Define list of file names (or dynamically fetch file list if applicable)
file_names = {
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CABLE-POP_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLASSIC_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLM6.0_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\EDv3_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\E3SM_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\IBIS_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ISBA-CTRIP_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\JSBACH_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPJmL_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPX-Bern_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\OCN_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ORCHIDEE_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\SDGVM_S3_ra.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\VISIT_S3_ra.nc",  
    };

% 定义时间范围
years = 2015:2023; % 提取的年份范围
num_years = length(years);
num_months = 12;

% 动态读取第一个文件的维度
ra_sample = ncread(file_names{1}, 'ra');
[lon_size, lat_size, ~] = size(ra_sample);

% 初始化矩阵，存储月数据和年数据
annual_ra_combined = nan(180, 360, num_months, num_years, length(file_names)); % 月数据
annual_ra_sum = nan(180, 360, num_years, length(file_names));                 % 年数据

% 定义每个月的天数（平年为基准）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% 遍历每个文件
for f = 1:length(file_names)
    f
    % 读取变量 'ra'
    ra_data = ncread(file_names{f}, 'ra'); % 假设 ra 数据是 3D 的：lon x lat x time
    [~, ~, time_dim] = size(ra_data);

    % 特殊处理 ISBA 文件
    if contains(file_names{f}, 'ISBA')
        % ISBA 文件：360x150 -> 150x360，添加 NaN 补齐为 180x360
        ra_data_rotated = rot90(ra_data, 1); % 旋转为 150x360
        ra_data_rotated = cat(1, ra_data_rotated, nan(30, 360, time_dim)); % 补齐为 180x360
    elseif contains(file_names{f}, 'CLM6')
        ra_data_rotated = rot90(ra_data, 1);
        ra_data_rotated=[ra_data_rotated(:,181:360,:),ra_data_rotated(:,1:180,:)];
    else
        % 普通文件：360x180 -> 180x360
        ra_data_rotated = rot90(ra_data, 1);
    end

    % 确定提取的时间范围
    start_idx = time_dim - num_years * num_months + 1; % 从 2015 年开始
    end_idx = time_dim;                               % 到 2023 年最后一个月

    % 提取 2015-2023 年的数据
    ra_sub_data = ra_data_rotated(:, :, start_idx:end_idx); % 提取最后 9 年的数据

    % 处理数据（按月累计并考虑闰年）
    for y = 1:num_years
        current_year = years(y);

        % 检查是否为闰年
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % 闰年 2 月有 29 天
        else
            days_in_month(2) = 28; % 平年 2 月有 28 天
        end

        % 计算当年的每月秒数
        seconds_per_month = days_in_month * 24 * 3600;

        % 提取当年数据
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        ra_yearly_data = ra_sub_data(:, :, month_start_idx:month_end_idx);

        % 转换单位并累计
        ra_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            ra_monthly(:, :, m) = ra_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % 存入月度矩阵
        annual_ra_combined(:, :, :, y, f) = ra_monthly;

        % 将 12 个月的数据求和，存入年度矩阵
        annual_ra_sum(:, :, y, f) = sum(ra_monthly, 3); % 对第 3 维求和
    end
end

%% TRENDY Rh

% Define list of file names (or dynamically fetch file list if applicable)
file_names = {
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CABLE-POP_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLASSIC_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLM6.0_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\EDv3_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\E3SM_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\IBIS_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ISBA-CTRIP_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\JSBACH_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPJmL_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPX-Bern_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\OCN_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ORCHIDEE_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\SDGVM_S3_rh.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\VISIT_S3_rh.nc",  
    };

% 定义时间范围
years = 2015:2023; % 提取的年份范围
num_years = length(years);
num_months = 12;

% 动态读取第一个文件的维度
rh_sample = ncread(file_names{1}, 'rh');
[lon_size, lat_size, ~] = size(rh_sample);

% 初始化矩阵，存储月数据和年数据
annual_rh_combined = nan(180, 360, num_months, num_years, length(file_names)); % 月数据
annual_rh_sum = nan(180, 360, num_years, length(file_names));                 % 年数据

% 定义每个月的天数（平年为基准）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% 遍历每个文件
for f = 1:length(file_names)
    f
    % 读取变量 'rh'
    rh_data = ncread(file_names{f}, 'rh'); % 假设 rh 数据是 3D 的：lon x lat x time
    [~, ~, time_dim] = size(rh_data);

    % 特殊处理 ISBA 文件
    if contains(file_names{f}, 'ISBA')
        % ISBA 文件：360x150 -> 150x360，添加 NaN 补齐为 180x360
        rh_data_rotated = rot90(rh_data, 1); % 旋转为 150x360
        rh_data_rotated = cat(1, rh_data_rotated, nan(30, 360, time_dim)); % 补齐为 180x360
    elseif contains(file_names{f}, 'CLM6')
        rh_data_rotated = rot90(rh_data, 1);
        rh_data_rotated=[rh_data_rotated(:,181:360,:),rh_data_rotated(:,1:180,:)];
    else
        % 普通文件：360x180 -> 180x360
        rh_data_rotated = rot90(rh_data, 1);
    end

    % 确定提取的时间范围
    start_idx = time_dim - num_years * num_months + 1; % 从 2015 年开始
    end_idx = time_dim;                               % 到 2023 年最后一个月

    % 提取 2015-2023 年的数据
    rh_sub_data = rh_data_rotated(:, :, start_idx:end_idx); % 提取最后 9 年的数据

    % 处理数据（按月累计并考虑闰年）
    for y = 1:num_years
        current_year = years(y);

        % 检查是否为闰年
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % 闰年 2 月有 29 天
        else
            days_in_month(2) = 28; % 平年 2 月有 28 天
        end

        % 计算当年的每月秒数
        seconds_per_month = days_in_month * 24 * 3600;

        % 提取当年数据
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        rh_yearly_data = rh_sub_data(:, :, month_start_idx:month_end_idx);

        % 转换单位并累计
        rh_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            rh_monthly(:, :, m) = rh_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % 存入月度矩阵
        annual_rh_combined(:, :, :, y, f) = rh_monthly;

        % 将 12 个月的数据求和，存入年度矩阵
        annual_rh_sum(:, :, y, f) = sum(rh_monthly, 3); % 对第 3 维求和
    end
end
%% TRENDY tas

% Define list of file names (or dynamically fetch file list if applicable)
file_names = {
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CABLE-POP_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLASSIC_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLM6.0_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\EDv3_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\E3SM_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\IBIS_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ISBA-CTRIP_S3_tasLut.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\JSBACH_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPJmL_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPX-Bern_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\OCN_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ORCHIDEE_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\SDGVM_S3_tas.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\VISIT_S3_tas.nc",
    };

% 定义时间范围
years = 2015:2023; % 提取的年份范围
num_years = length(years);
num_months = 12;

% 动态读取第一个文件的维度
tas_sample = ncread(file_names{1}, 'tas');
[lon_size, lat_size, ~] = size(tas_sample);

% 初始化矩阵，存储月数据和年数据
annual_tas_combined = nan(180, 360, num_months, num_years, length(file_names)); % 月数据
annual_tas_sum = nan(180, 360, num_years, length(file_names));                 % 年数据

% 定义每个月的天数（平年为基准）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% 遍历每个文件
for f = 1:length(file_names)
    f
    % 读取变量 'tas'
    if contains(file_names{f}, 'ISBA')
        tas_data = ncread(file_names{f}, 'tasLut'); % 假设 tas 数据是 3D 的：lon x lat x time
        [~, ~, time_dim] = size(tas_data);
    else
        tas_data = ncread(file_names{f}, 'tas'); % 假设 tas 数据是 3D 的：lon x lat x time
        [~, ~, time_dim] = size(tas_data);
    end

    % 特殊处理 ISBA 文件
    if contains(file_names{f}, 'ISBA')
        % ISBA 文件：360x150 -> 150x360，添加 NaN 补齐为 180x360
        tas_data_rotated = rot90(tas_data, 1); % 旋转为 150x360
        tas_data_rotated = cat(1, tas_data_rotated, nan(30, 360, time_dim)); % 补齐为 180x360
    elseif contains(file_names{f}, 'CLM6')
        tas_data_rotated = rot90(tas_data, 1);
        tas_data_rotated=[tas_data_rotated(:,181:360,:),tas_data_rotated(:,1:180,:)];
    else
        % 普通文件：360x180 -> 180x360
        tas_data_rotated = rot90(tas_data, 1);
    end

    % 确定提取的时间范围
    start_idx = time_dim - num_years * num_months + 1; % 从 2015 年开始
    end_idx = time_dim;                               % 到 2023 年最后一个月

    % 提取 2015-2023 年的数据
    tas_sub_data = tas_data_rotated(:, :, start_idx:end_idx); % 提取最后 9 年的数据

    % 处理数据（按月累计并考虑闰年）
    for y = 1:num_years
        current_year = years(y);

        % 提取当年数据
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        tas_yearly_data = tas_sub_data(:, :, month_start_idx:month_end_idx);

        % 转换单位并累计
        tas_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            tas_monthly(:, :, m) = tas_yearly_data(:, :, m)-273.15;
        end

        % 存入月度矩阵
        annual_tas_combined(:, :, :, y, f) = tas_monthly;

        % 将 12 个月的数据求和，存入年度矩阵
        annual_tas_sum(:, :, y, f) = mean(tas_monthly, 3); % 对第 3 维求和
    end
end
%% TRENDY mrso

% Define list of file names (or dynamically fetch file list if applicable)
file_names = {
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CABLE-POP_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLASSIC_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\CLM6.0_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\EDv3_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\E3SM_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\IBIS_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ISBA-CTRIP_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\JSBACH_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPJmL_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\LPX-Bern_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\OCN_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\ORCHIDEE_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\SDGVM_S3_mrso.nc",
    "E:\phd_file\Boreal_North_America\TRENDY\aggregate\VISIT_S3_mrso.nc",
    };

% 定义时间范围
years = 2015:2023; % 提取的年份范围
num_years = length(years);
num_months = 12;

% 动态读取第一个文件的维度
mrso_sample = ncread(file_names{1}, 'mrso');
[lon_size, lat_size, ~] = size(mrso_sample);

% 初始化矩阵，存储月数据和年数据
annual_mrso_combined = nan(180, 360, num_months, num_years, length(file_names)); % 月数据
annual_mrso_sum = nan(180, 360, num_years, length(file_names));                 % 年数据

% 定义每个月的天数（平年为基准）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% 遍历每个文件
for f = 1:length(file_names)
    f
    % 读取变量 'mrso'
    mrso_data = ncread(file_names{f}, 'mrso'); % 假设 mrso 数据是 3D 的：lon x lat x time
    [~, ~, time_dim] = size(mrso_data);


    % 特殊处理 ISBA 文件
    if contains(file_names{f}, 'ISBA')
        % ISBA 文件：360x150 -> 150x360，添加 NaN 补齐为 180x360
        mrso_data_rotated = rot90(mrso_data, 1); % 旋转为 150x360
        mrso_data_rotated = cat(1, mrso_data_rotated, nan(30, 360, time_dim)); % 补齐为 180x360
    elseif contains(file_names{f}, 'CLM6')
        mrso_data_rotated = rot90(mrso_data, 1);
        mrso_data_rotated=[mrso_data_rotated(:,181:360,:),mrso_data_rotated(:,1:180,:)];
    else
        % 普通文件：360x180 -> 180x360
        mrso_data_rotated = rot90(mrso_data, 1);
    end

    % 确定提取的时间范围
    start_idx = time_dim - num_years * num_months + 1; % 从 2015 年开始
    end_idx = time_dim;                               % 到 2023 年最后一个月

    % 提取 2015-2023 年的数据
    mrso_sub_data = mrso_data_rotated(:, :, start_idx:end_idx); % 提取最后 9 年的数据

    % 处理数据（按月累计并考虑闰年）
    for y = 1:num_years
        current_year = years(y);

        % 提取当年数据
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        mrso_yearly_data = mrso_sub_data(:, :, month_start_idx:month_end_idx);

        % 转换单位并累计
        mrso_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            mrso_monthly(:, :, m) = mrso_yearly_data(:, :, m);
        end

        % 存入月度矩阵
        annual_mrso_combined(:, :, :, y, f) = mrso_monthly;

        % 将 12 个月的数据求和，存入年度矩阵
        annual_mrso_sum(:, :, y, f) = mean(mrso_monthly, 3); % 对第 3 维求和
    end
end

save TRENDY_variable.mat annual_GPP_combined annual_GPP_sum annual_mrso_combined annual_mrso_sum annual_ra_combined annual_ra_sum annual_rh_combined ...
    annual_rh_sum annual_tas_combined annual_tas_sum

