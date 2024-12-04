%% Mask create
clc
clear
% Canada mask
pixel_mask=importdata('E:\phd_file\Boreal_North_America\Canada_mask.tif');
forest_mask=importdata("E:\phd_file\Boreal_North_America\region_lu.tif");
forest_mask(forest_mask~=1)=nan;
pixel_mask=pixel_mask.*forest_mask;

% 创建环境
data=geotiffread('E:\phd_file\Boreal_North_America\region_lu.tif');
info=geotiffinfo('E:\phd_file\Boreal_North_America\region_lu.tif');
[m,n] = size(data);
k=1;
for i = 1:m
    for j = 1:1
        [lat,lon]= pix2latlon(info.RefMatrix, i, j);   %读取栅格数据第1列所有行的纬度；
        lat_(k,:)=lat; %将纬度数据存储为1列；
        k=k+1;
    end
end

k=1;
for ii = 1:1
    for jj = 1:n
        [lat,lon]= pix2latlon(info.RefMatrix, ii, jj);   %读取栅格数据第1行所有行的经度；
        lon_(k,:)=lon;  %将经度数据存储为1列；
        k=k+1;
    end
end
[lon1,lat1]=meshgrid(lon_,lat_);

Boundry = shaperead("E:\phd_file\yuling_shiliang\countries.shp");
bou_canX = [Boundry(:).X];
bou_canY= [Boundry(:).Y];

clearvars -except lon1 lat1 pixel_mask bou_canX bou_canY m n
%% CMIP6 GPP

% Define list of file names (or dynamically fetch file list if applicable)
GPP_file_names = {
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_BCC-CSM2-MR_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_CESM2-WACCM_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_CMCC-CM2-SR5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_GFDL-ESM4_ssp126_r1i1p1f1_gr1_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gr_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_MPI-ESM1-2-LR_ssp126_r1i1p1f1_gn_201501-203412.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_NorESM2-LM_ssp126_r1i1p1f1_gn_202101-203012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_NorESM2-MM_ssp126_r1i1p1f1_gn_201502-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_NorESM2-MM_ssp126_r1i1p1f1_gn_202101-203012.nc",  
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\gpp_Lmon_TaiESM1_ssp126_r1i1p1f1_gn_201501-210012.nc",
    };

% Define grid dimensions and time variables
years = 2015:2030;
summer_months = [6, 7, 8]; % June, July, August
seconds_per_month = [30*24*3600, 31*24*3600, 31*24*3600]; % seconds in June, July, August


% Initialize 4D matrix to store results for all files
summer_GPP_all = nan(180, 360, length(years), length(GPP_file_names)-2);


% Loop over each year
f = 1; % Initialize f outside the loop
k=1;
while f <= length(GPP_file_names)
% f
    if contains(GPP_file_names{f}, 'NorESM2-MM') && contains(GPP_file_names{f}, '201502-202012')

        % Read data from 2015 February to 2020 December file
        GPP_data_part1 = ncread(GPP_file_names{f}, 'gpp');
        GPP_data_part2 = ncread(GPP_file_names{f+1}, 'gpp');

        % Concatenate data across the time dimension for a full series
        GPP_data = cat(3, GPP_data_part1, GPP_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = -1;

    elseif contains(GPP_file_names{f}, 'NorESM2-LM') && contains(GPP_file_names{f}, '201501-202012')

        % Read data from 2015 February to 2020 December file
        GPP_data_part1 = ncread(GPP_file_names{f}, 'gpp');
        GPP_data_part2 = ncread(GPP_file_names{f+1}, 'gpp');

        % Concatenate data across the time dimension for a full series
        GPP_data = cat(3, GPP_data_part1, GPP_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = 0;

    else
        % For other models, read normally (January 2015 start)

        GPP_data = ncread(GPP_file_names{f}, 'gpp');
        start_offset = 0;
    end

    % Initialize matrix to store summer GPP for each year
    summer_GPP = nan(360, 180, length(years));

    % Loop through each year to extract and process summer data
    for i = 1:length(years)
        % Calculate starting index with offset adjustment for February start
        start_index = (years(i) - 2015) * 12 + summer_months(1) + start_offset;
        indices = start_index:start_index + 2;

        % Adjust each month's GPP by seconds in that month
        monthly_GPP_adjusted = GPP_data(:, :, indices) .* reshape(seconds_per_month, [1, 1, length(seconds_per_month)]);

        % Sum the adjusted GPP data over the summer months
        summer_GPP(:, :, i) = sum(monthly_GPP_adjusted, 3);
    end

    % Convert to g/m^2 by multiplying by 1000
    summer_GPP = summer_GPP * 1000;

    % Rotate and re-arrange each layer for final dimensions
    rotated_GPP = nan(180, 360, length(years));
    for i = 1:size(summer_GPP, 3)
        % Rotate 90° counterclockwise and adjust longitude order
        GPP_temp = rot90(summer_GPP(:, :, i), 1);
        GPP_temp = [GPP_temp(:, 181:360), GPP_temp(:, 1:180)];
        rotated_GPP(:, :, i) = GPP_temp.*pixel_mask;
    end

    % Store in the 4D matrix
    summer_GPP_all(:, :, :, k) = rotated_GPP;
    k=k+1;
    % Move to the next file
    f = f + 1;
end
%% CMIP6 NEP
% Define list of file names (or dynamically fetch file list if applicable)
nep_file_names = {
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_BCC-CSM2-MR_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_CESM2-WACCM_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_CMCC-CM2-SR5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_GFDL-ESM4_ssp126_r1i1p1f1_gr1_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gr_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_MPI-ESM1-2-LR_ssp126_r1i1p1f1_gn_201501-203412.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_NorESM2-LM_ssp126_r1i1p1f1_gn_202101-203012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_NorESM2-MM_ssp126_r1i1p1f1_gn_201502-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_NorESM2-MM_ssp126_r1i1p1f1_gn_202101-203012.nc",      
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\nep_Emon_TaiESM1_ssp126_r1i1p1f1_gn_201501-210012.nc"};

% Define grid dimensions and time variables
years = 2015:2030;
summer_months = [6, 7, 8]; % June, July, August
seconds_per_month = [30*24*3600, 31*24*3600, 31*24*3600]; % seconds in June, July, August


% Initialize 4D matrix to store results for all files
summer_nep_all = nan(180, 360, length(years), length(nep_file_names)-2);


% Loop over each year
f = 1; % Initialize f outside the loop
k=1;
while f <= length(nep_file_names)

    if contains(nep_file_names{f}, 'NorESM2-MM') && contains(nep_file_names{f}, '201502-202012')

        % Read data from 2015 February to 2020 December file
        nep_data_part1 = ncread(nep_file_names{f}, 'nep');
        nep_data_part2 = ncread(nep_file_names{f+1}, 'nep');

        % Concatenate data across the time dimension for a full series
        nep_data = cat(3, nep_data_part1, nep_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = -1;

    elseif contains(nep_file_names{f}, 'NorESM2-LM') && contains(nep_file_names{f}, '201501-202012')

        % Read data from 2015 February to 2020 December file
        nep_data_part1 = ncread(nep_file_names{f}, 'nep');
        nep_data_part2 = ncread(nep_file_names{f+1}, 'nep');

        % Concatenate data across the time dimension for a full series
        nep_data = cat(3, nep_data_part1, nep_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = 0;

    else
        % For other models, read normally (January 2015 start)

        nep_data = ncread(nep_file_names{f}, 'nep');
        start_offset = 0;
    end

    % Initialize matrix to store summer nep for each year
    summer_nep = nan(360, 180, length(years));

    % Loop through each year to extract and process summer data
    for i = 1:length(years)
        % Calculate starting index with offset adjustment for February start
        start_index = (years(i) - 2015) * 12 + summer_months(1) + start_offset;
        indices = start_index:start_index + 2;

        % Adjust each month's nep by seconds in that month
        monthly_nep_adjusted = nep_data(:, :, indices) .* reshape(seconds_per_month, [1, 1, length(seconds_per_month)]);

        % Sum the adjusted nep data over the summer months
        summer_nep(:, :, i) = sum(monthly_nep_adjusted, 3);
    end

    % Convert to g/m^2 by multiplying by 1000
    summer_nep = summer_nep * 1000;

    % Rotate and re-arrange each layer for final dimensions
    rotated_nep = nan(180, 360, length(years));
    for i = 1:size(summer_nep, 3)
        % Rotate 90° counterclockwise and adjust longitude order
        nep_temp = rot90(summer_nep(:, :, i), 1);
        nep_temp = [nep_temp(:, 181:360), nep_temp(:, 1:180)];
        rotated_nep(:, :, i) = nep_temp.*pixel_mask;
    end

    % Store in the 4D matrix
    summer_nep_all(:, :, :, k) = rotated_nep;
    k=k+1;
    % Move to the next file
    f = f + 1;
end
%% CMIP6 ER
for i=1:size(summer_nep_all,4)
    for j=1:size(summer_nep_all,3)
        nep_temp=summer_nep_all(:, :, j, i);
        gpp_temp=summer_GPP_all(:, :, j, i);
        summer_ER_all(:, :, j, i)=gpp_temp-nep_temp;
    end
end
%% CMIP6 SM
% Define list of file names (or dynamically fetch file list if applicable)
mrso_file_names = {
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_BCC-CSM2-MR_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_CESM2-WACCM_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_CMCC-CM2-SR5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_GFDL-ESM4_ssp126_r1i1p1f1_gr1_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gr_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_MPI-ESM1-2-LR_ssp126_r1i1p1f1_gn_201501-203412.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_NorESM2-LM_ssp126_r1i1p1f1_gn_202101-203012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_NorESM2-MM_ssp126_r1i1p1f1_gn_201502-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_NorESM2-MM_ssp126_r1i1p1f1_gn_202101-203012.nc",    
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\mrso_Lmon_TaiESM1_ssp126_r1i1p1f1_gn_201501-210012.nc"};

% Define grid dimensions and time variables
years = 2015:2030;
summer_months = [6, 7, 8]; % June, July, August
% seconds_per_month = [30*24*3600, 31*24*3600, 31*24*3600]; % seconds in June, July, August


% Initialize 4D matrix to store results for all files
summer_mrso_all = nan(180, 360, length(years), length(mrso_file_names)-2);


% Loop over each year
f = 1; % Initialize f outside the loop
k=1;
while f <= length(mrso_file_names)

    if contains(mrso_file_names{f}, 'NorESM2-MM') && contains(mrso_file_names{f}, '201502-202012')

        % Read data from 2015 February to 2020 December file
        mrso_data_part1 = ncread(mrso_file_names{f}, 'mrso');
        mrso_data_part2 = ncread(mrso_file_names{f+1}, 'mrso');

        % Concatenate data across the time dimension for a full series
        mrso_data = cat(3, mrso_data_part1, mrso_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = -1;

    elseif contains(mrso_file_names{f}, 'NorESM2-LM') && contains(mrso_file_names{f}, '201501-202012')

        % Read data from 2015 February to 2020 December file
        mrso_data_part1 = ncread(mrso_file_names{f}, 'mrso');
        mrso_data_part2 = ncread(mrso_file_names{f+1}, 'mrso');

        % Concatenate data across the time dimension for a full series
        mrso_data = cat(3, mrso_data_part1, mrso_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = 0;

    else
        % For other models, read normally (January 2015 start)

        mrso_data = ncread(mrso_file_names{f}, 'mrso');
        start_offset = 0;
    end

    % Initialize matrix to store summer mrso for each year
    summer_mrso = nan(360, 180, length(years));

    % Loop through each year to extract and process summer data
    for i = 1:length(years)
        % Calculate starting index with offset adjustment for February start
        start_index = (years(i) - 2015) * 12 + summer_months(1) + start_offset;
        indices = start_index:start_index + 2;

        % Adjust each month's mrso by seconds in that month
        monthly_mrso_adjusted = mrso_data(:, :, indices);

        % Sum the adjusted mrso data over the summer months
        summer_mrso(:, :, i) = mean(monthly_mrso_adjusted, 3);
    end


    % Rotate and re-arrange each layer for final dimensions
    rotated_mrso = nan(180, 360, length(years));
    for i = 1:size(summer_mrso, 3)
        % Rotate 90° counterclockwise and adjust longitude order
        mrso_temp = rot90(summer_mrso(:, :, i), 1);
        mrso_temp = [mrso_temp(:, 181:360), mrso_temp(:, 1:180)];
        rotated_mrso(:, :, i) = mrso_temp.*pixel_mask;
    end

    % Store in the 4D matrix
    summer_mrso_all(:, :, :, k) = rotated_mrso;
    k=k+1;
    % Move to the next file
    f = f + 1;
end
%% CMIP6 TEM
% Define list of file names (or dynamically fetch file list if applicable)
tas_file_names = {
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_BCC-CSM2-MR_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_CanESM5_ssp126_r1i1p1f1_gn_201501-210012.nc"
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_CESM2-WACCM_ssp126_r1i1p1f1_gn_201501-206412.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_CMCC-CM2-SR5_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_GFDL-ESM4_ssp126_r1i1p1f1_gr1_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_IPSL-CM6A-LR_ssp126_r1i1p1f1_gr_201501-210012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_MPI-ESM1-2-LR_ssp126_r1i1p1f1_gn_201501-203412.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_NorESM2-LM_ssp126_r1i1p1f1_gn_202101-203012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_NorESM2-MM_ssp126_r1i1p1f1_gn_201501-202012.nc",
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_NorESM2-MM_ssp126_r1i1p1f1_gn_202101-203012.nc",     
    "E:\phd_file\Boreal_North_America\CMIP6\aggregate\tas_Amon_TaiESM1_ssp126_r1i1p1f1_gn_201501-210012.nc"};

% Define grid dimensions and time variables
years = 2015:2030;
summer_months = [6, 7, 8]; % June, July, August
seconds_per_month = [30*24*3600, 31*24*3600, 31*24*3600]; % seconds in June, July, August


% Initialize 4D matrix to store results for all files
summer_tas_all = nan(180, 360, length(years), length(tas_file_names)-2);


% Loop over each year
f = 1; % Initialize f outside the loop
k=1;
while f <= length(tas_file_names)

    if contains(tas_file_names{f}, 'NorESM2-MM') && contains(tas_file_names{f}, '201501-202012')

        % Read data from 2015 February to 2020 December file
        tas_data_part1 = ncread(tas_file_names{f}, 'tas');
        tas_data_part2 = ncread(tas_file_names{f+1}, 'tas');

        % Concatenate data across the time dimension for a full series
        tas_data = cat(3, tas_data_part1, tas_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = 0;

    elseif contains(tas_file_names{f}, 'NorESM2-LM') && contains(tas_file_names{f}, '201501-202012')

        % Read data from 2015 February to 2020 December file
        tas_data_part1 = ncread(tas_file_names{f}, 'tas');
        tas_data_part2 = ncread(tas_file_names{f+1}, 'tas');

        % Concatenate data across the time dimension for a full series
        tas_data = cat(3, tas_data_part1, tas_data_part2);
        f = f + 1; % Skip the next file in loop since it's already used here

        % Adjust start month to February 2015 by using an offset of -1
        start_offset = 0;

    else
        % For other models, read normally (January 2015 start)

        tas_data = ncread(tas_file_names{f}, 'tas');
        start_offset = 0;
    end

    % Initialize matrix to store summer tas for each year
    summer_tas = nan(360, 180, length(years));

    % Loop through each year to extract and process summer data
    for i = 1:length(years)
        % Calculate starting index with offset adjustment for February start
        start_index = (years(i) - 2015) * 12 + summer_months(1) + start_offset;
        indices = start_index:start_index + 2;

        % Adjust each month's tas by seconds in that month
        monthly_tas_adjusted = tas_data(:, :, indices);

        % Sum the adjusted tas data over the summer months
        summer_tas(:, :, i) = mean(monthly_tas_adjusted, 3)-273.15;
    end



    % Rotate and re-arrange each layer for final dimensions
    rotated_tas = nan(180, 360, length(years));
    for i = 1:size(summer_tas, 3)
        % Rotate 90° counterclockwise and adjust longitude order
        tas_temp = rot90(summer_tas(:, :, i), 1);
        tas_temp = [tas_temp(:, 181:360), tas_temp(:, 1:180)];
        rotated_tas(:, :, i) = tas_temp.*pixel_mask;
    end

    % Store in the 4D matrix
    summer_tas_all(:, :, :, k) = rotated_tas;
    k=k+1;
    % Move to the next file
    f = f + 1;
end
%%

for i=1:size(summer_ER_all,4)
    for j=1:size(summer_ER_all,3)


        ER_temp=summer_ER_all(:, :, j, i);
        GPP_temp=summer_GPP_all(:, :, j, i);
        NEP_temp=summer_nep_all(:, :, j, i);
        TEM_temp=summer_tas_all(:, :, j, i);
        SM_temp=summer_mrso_all(:, :, j, i);

        ER_temp(ER_temp==0)=nan;
        GPP_temp(GPP_temp==0)=nan;
        NEP_temp(NEP_temp==0)=nan;
        TEM_temp(TEM_temp==0)=nan;
        SM_temp(SM_temp==0)=nan;
        area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
        area_grid(SM_temp==0)=nan;

        ER_list(i,j)=nansum(nansum(ER_temp.*area_grid/(10^15)));
        GPP_list(i,j)=nansum(nansum(GPP_temp.*area_grid/(10^15)));
        NEP_list(i,j)=nansum(nansum(NEP_temp.*area_grid/(10^15)));
        TEM_list(i,j)=nansum(nansum(TEM_temp.*area_grid))/(nansum(nansum(area_grid)));
        SM_list(i,j)=nansum(nansum(SM_temp.*area_grid))/(nansum(nansum(area_grid)));

    end
end
%%
filelist={"ACCESS-ESM1-5","BCC-CSM2-MR","CanESM5","CESM2-WACCM","CMCC-CM2-SR5","CMCC-ESM2","GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM","TaiESM1",};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
ii=1;
for i=1:6
    nexttile
    scatter(SM_list(i,:),ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    % 使用二项式拟合 (二次多项式)
    [p,s] = polyfit(SM_list(i,:),ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(SM_list(i,:)),max(SM_list(i,:)));
    % 计算拟合值
    [y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    % 绘制置信区间阴影
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['R^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.074484185627727 0],'FontName','Arial','FontSize',12,'Color','k')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('RZSM (kg m^{-2})','FontName','Arial','FontSize',12)
    ylabel('ER (PgC)','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',12,'fontweight','bold')

    nexttile
    scatter(TEM_list(i,:),ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(TEM_list(i,:),ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(TEM_list(i,:)),max(TEM_list(i,:)));
    % 计算拟合值
    [y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    % 绘制置信区间阴影
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['R^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.074484185627727 0],'FontName','Arial','FontSize',12,'Color','k')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('TEM (°C)','FontName','Arial','FontSize',12)
    ylabel('ER (PgC)','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii+1},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',12,'fontweight','bold')
    ii=ii+2;
end

% Assuming your subplots have already been created in a figure
figure_handle = gcf;  % Get current figure handle
axes_handles = findall(figure_handle, 'type', 'axes');  % Get all axes handles

% Define a padding factor to add extra space around scatter points
padding_factor = 0.1;  % Adjust padding as needed

for i = 1:length(axes_handles)
    ax = axes_handles(i);  % Get current axis
    scatter_data = findobj(ax, 'Type', 'Scatter');  % Find scatter plot in the current axis

    if ~isempty(scatter_data)
        % Get x and y data from scatter plot
        xData = scatter_data.XData;
        yData = scatter_data.YData;

        % Calculate axis limits with padding
        xRange = max(xData) - min(xData);
        yRange = max(yData) - min(yData);

        xlim(ax, [min(xData) - padding_factor * xRange, max(xData) + padding_factor * xRange]);
        ylim(ax, [min(yData) - padding_factor * yRange, max(yData) + padding_factor * yRange]);
    end
end

set(gcf,'unit','centimeters','position',[12.832291666666668,13.599583333333335,30.718125,21.246041666666677]);
result=['E:\phd_file\Boreal_North_America\Result\V6\CMIP6_indiviual_scatter1.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"ACCESS-ESM1-5","BCC-CSM2-MR","CanESM5","CESM2-WACCM","CMCC-CM2-SR5","CMCC-ESM2","GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM","TaiESM1",};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
ii=1;
for i=7:12
    nexttile
    scatter(SM_list(i,:),ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    % 使用二项式拟合 (二次多项式)
    [p,s] = polyfit(SM_list(i,:),ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(SM_list(i,:)),max(SM_list(i,:)));
    % 计算拟合值
    [y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    % 绘制置信区间阴影
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['R^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.074484185627727 0],'FontName','Arial','FontSize',12,'Color','k')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('RZSM (kg m^{-2})','FontName','Arial','FontSize',12)
    ylabel('ER (PgC)','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',12,'fontweight','bold')

    nexttile
    scatter(TEM_list(i,:),ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(TEM_list(i,:),ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(TEM_list(i,:)),max(TEM_list(i,:)));
    % 计算拟合值
    [y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    % 绘制置信区间阴影
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['R^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.074484185627727 0],'FontName','Arial','FontSize',12,'Color','k')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('TEM (°C)','FontName','Arial','FontSize',12)
    ylabel('ER (PgC)','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii+1},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',12,'fontweight','bold')
    ii=ii+2;
end

% Assuming your subplots have already been created in a figure
figure_handle = gcf;  % Get current figure handle
axes_handles = findall(figure_handle, 'type', 'axes');  % Get all axes handles

% Define a padding factor to add extra space around scatter points
padding_factor = 0.1;  % Adjust padding as needed

for i = 1:length(axes_handles)
    ax = axes_handles(i);  % Get current axis
    scatter_data = findobj(ax, 'Type', 'Scatter');  % Find scatter plot in the current axis

    if ~isempty(scatter_data)
        % Get x and y data from scatter plot
        xData = scatter_data.XData;
        yData = scatter_data.YData;

        % Calculate axis limits with padding
        xRange = max(xData) - min(xData);
        yRange = max(yData) - min(yData);

        xlim(ax, [min(xData) - padding_factor * xRange, max(xData) + padding_factor * xRange]);
        ylim(ax, [min(yData) - padding_factor * yRange, max(yData) + padding_factor * yRange]);
    end
end

set(gcf,'unit','centimeters','position',[12.832291666666668,13.599583333333335,30.718125,21.246041666666677]);
result=['E:\phd_file\Boreal_North_America\Result\V6\CMIP6_indiviual_scatter2.png']
% print(result,ff,'-r600','-dpng' );
%% 计算后把结果放入TRENDY_scatter中出图
summer_ER_all=nanmean(summer_ER_all,4);
summer_GPP_all=nanmean(summer_GPP_all,4);
summer_nep_all=nanmean(summer_nep_all,4);
summer_tas_all=nanmean(summer_tas_all,4);
summer_mrso_all=nanmean(summer_mrso_all,4);


for j=1:size(summer_ER_all,3)


    ER_temp=summer_ER_all(:, :, j);
    GPP_temp=summer_GPP_all(:, :, j);
    NEP_temp=summer_nep_all(:, :, j);
    TEM_temp=summer_tas_all(:, :, j);
    SM_temp=summer_mrso_all(:, :, j);

    ER_temp(ER_temp==0)=nan;
    GPP_temp(GPP_temp==0)=nan;
    NEP_temp(NEP_temp==0)=nan;
    TEM_temp(TEM_temp==0)=nan;
    SM_temp(SM_temp==0)=nan;
    area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
    area_grid(SM_temp==0)=nan;

    ER_list2(j)=nansum(nansum(ER_temp.*area_grid/(10^15)));
    GPP_list2(j)=nansum(nansum(GPP_temp.*area_grid/(10^15)));
    NEP_list2(j)=nansum(nansum(NEP_temp.*area_grid/(10^15)));
    TEM_list2(j)=nansum(nansum(TEM_temp.*area_grid))/(nansum(nansum(area_grid)));
    SM_list2(j)=nansum(nansum(SM_temp.*area_grid))/(nansum(nansum(area_grid)));

end
% save CMIP6_summer_variable_result.mat ER_list2 GPP_list2 NEP_list2 TEM_list2 SM_list2
