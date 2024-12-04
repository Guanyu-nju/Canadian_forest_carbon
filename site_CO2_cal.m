clc
clear
time=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_etl_surface-insitu_6_allvalid.nc",'time_components');
% 转换为datetime对象，注意矩阵中的每一列是一个时间
utc_times = datetime(time(1,:), time(2,:), time(3,:), time(4,:), time(5,:), time(6,:), 'TimeZone', 'UTC');
% 将时间转换为UTC-6 (加拿大当地时间)
local_times = utc_times;
local_times.TimeZone = 'America/Winnipeg';  % 例如使用UTC-6的埃德蒙顿时区
daytime_mask=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_etl_surface-insitu_6_allvalid.nc",'CT_assim');
daytime_mask=double(daytime_mask);
daytime_mask(daytime_mask~=1)=nan;
value=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_etl_surface-insitu_6_allvalid.nc",'value')*10^6;
value(value>450)=nan;
value = sgolayfilt(value,2,21);     % sg滤波，参数设置为默认的3,5
daytime_value=value.*daytime_mask;

% 筛选出时间范围内的数据
mask = (year(local_times) >= 2015) & (year(local_times) <= 2023);
filtered_times = local_times(mask)';
filtered_concentration = daytime_value(mask);

% 创建年和月列，作为分组变量
years = year(filtered_times);   % 提取年份
months = month(filtered_times); % 提取月份
% 按年和月分组，并计算每个月的平均浓度
% 将数据转换为表格形式
data_table = table(filtered_times, filtered_concentration, years, months);

% 使用 groupsummary 按年和月份分组并计算平均浓度
monthly_avg = groupsummary(data_table, {'years','months'}, 'mean', 'filtered_concentration');


% 进一步筛选每天14:00到16:00的数据
timeMask = (hour(filtered_times) >= 0) & (hour(filtered_times) <= 4);
night_values=value(mask);
night_values=night_values(timeMask);
night_times=filtered_times(timeMask);
night_years = year(night_times);   % 提取年份
night_months = month(night_times); % 提取月份
night_data_table = table(night_times, night_values, night_years, night_months);
night_monthly_avg = groupsummary(night_data_table, {'night_years','night_months'}, 'mean', 'night_values');


% 将table输出为Excel文件
filename = 'E:\phd_file\Boreal_North_America\Ground_flux\Daytime_etl_monthly_co2_SG.xlsx'; % Excel文件名
writetable(monthly_avg, filename);
% filename = 'E:\phd_file\Boreal_North_America\Ground_flux\Nighttime_etl_monthly_co2.xlsx'; % Excel文件名
% writetable(night_monthly_avg, filename);

%%
% 显示每个月份的平均CO2浓度
cc=monthly_avg{:,4}
cc=detrend(cc)
k=1
for ii=1:9
    list(ii)=mean(cc(k:12*ii));

    ii
    k=k+12
end
%%
clc
clear
time=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_llb_surface-insitu_6_allvalid.nc",'time_components');
% 转换为datetime对象，注意矩阵中的每一列是一个时间
utc_times = datetime(time(1,:), time(2,:), time(3,:), time(4,:), time(5,:), time(6,:), 'TimeZone', 'UTC');
% 将时间转换为UTC-7 (加拿大当地时间)
local_times = utc_times;
local_times.TimeZone = 'America/Edmonton';  % 例如使用UTC-7的埃德蒙顿时区

daytime_mask=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_llb_surface-insitu_6_allvalid.nc",'CT_assim');
daytime_mask=double(daytime_mask);
daytime_mask(daytime_mask~=1)=nan;
value=ncread("E:\phd_file\Boreal_North_America\Ground_flux\co2_llb_surface-insitu_6_allvalid.nc",'value')*10^6;
value(value>510)=nan;
daytime_value=value.*daytime_mask;

% 筛选出时间范围内的数据
mask = (year(local_times) >= 2015) & (year(local_times) <= 2023);
filtered_times = local_times(mask)';
filtered_concentration = daytime_value(mask);

% 创建年和月列，作为分组变量
years = year(filtered_times);   % 提取年份
months = month(filtered_times); % 提取月份
% 按年和月分组，并计算每个月的平均浓度
% 将数据转换为表格形式
data_table = table(filtered_times, filtered_concentration, years, months);

% 使用 groupsummary 按年和月份分组并计算平均浓度
monthly_avg = groupsummary(data_table, {'years','months'}, 'mean', 'filtered_concentration');


% 进一步筛选每天14:00到16:00的数据
timeMask = (hour(filtered_times) >= 0) & (hour(filtered_times) <= 4);
night_values=value(mask);
night_values=night_values(timeMask);
night_times=filtered_times(timeMask);
night_years = year(night_times);   % 提取年份
night_months = month(night_times); % 提取月份
night_data_table = table(night_times, night_values, night_years, night_months);
night_monthly_avg = groupsummary(night_data_table, {'night_years','night_months'}, 'mean', 'night_values');


% % 将table输出为Excel文件
% filename = 'E:\phd_file\Boreal_North_America\Ground_flux\Daytime_llb_monthly_co2.xlsx'; % Excel文件名
% writetable(monthly_avg, filename);
% filename = 'E:\phd_file\Boreal_North_America\Ground_flux\Nighttime_llb_monthly_co2.xlsx'; % Excel文件名
% writetable(night_monthly_avg, filename);