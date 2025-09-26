clc
clear
%% TRENDY GPP Data Processing

% Define list of file names for GPP data
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

% Define time range for data extraction
years = 2015:2023; % Year range to extract
num_years = length(years);
num_months = 12;

% Dynamically read dimensions from first file
gpp_sample = ncread(file_names{1}, 'gpp');
[lon_size, lat_size, ~] = size(gpp_sample);

% Initialize matrices for monthly and annual data storage
annual_GPP_combined = nan(180, 360, num_months, num_years, length(file_names)); % Monthly data
annual_GPP_sum = nan(180, 360, num_years, length(file_names));                 % Annual data

% Define days in each month (non-leap year as baseline)
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Process each file
for f = 1:length(file_names)
    f % Display current file number
    % Read 'gpp' variable
    gpp_data = ncread(file_names{f}, 'gpp'); % Assume gpp data is 3D: lon x lat x time
    [~, ~, time_dim] = size(gpp_data);

    % Special handling for ISBA files
    if contains(file_names{f}, 'ISBA')
        % ISBA files: 360x150 -> 150x360, add NaN to pad to 180x360
        gpp_data_rotated = rot90(gpp_data, 1); % Rotate to 150x360
        gpp_data_rotated = cat(1, gpp_data_rotated, nan(30, 360, time_dim)); % Pad to 180x360
    elseif contains(file_names{f}, 'CLM6')
        % CLM6 files: special rotation and reorganization
        gpp_data_rotated = rot90(gpp_data, 1);
        gpp_data_rotated = [gpp_data_rotated(:,181:360,:), gpp_data_rotated(:,1:180,:)];
    else
        % Regular files: 360x180 -> 180x360
        gpp_data_rotated = rot90(gpp_data, 1);
    end

    % Determine extraction time range
    start_idx = time_dim - num_years * num_months + 1; % Start from 2015
    end_idx = time_dim;                               % End at last month of 2023

    % Extract 2015-2023 data
    gpp_sub_data = gpp_data_rotated(:, :, start_idx:end_idx); % Extract last 9 years of data

    % Process data (monthly accumulation with leap year consideration)
    for y = 1:num_years
        current_year = years(y);

        % Check if leap year
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % Leap year: February has 29 days
        else
            days_in_month(2) = 28; % Non-leap year: February has 28 days
        end

        % Calculate seconds per month for current year
        seconds_per_month = days_in_month * 24 * 3600;

        % Extract current year's data
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        gpp_yearly_data = gpp_sub_data(:, :, month_start_idx:month_end_idx);

        % Convert units and accumulate
        gpp_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            % Convert from kg/m²/s to gC/m²/month
            gpp_monthly(:, :, m) = gpp_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % Store in monthly matrix
        annual_GPP_combined(:, :, :, y, f) = gpp_monthly;

        % Sum 12 months of data and store in annual matrix
        annual_GPP_sum(:, :, y, f) = sum(gpp_monthly, 3); % Sum along 3rd dimension
    end
end


%% TRENDY Ra (Autotrophic Respiration) Data Processing

% Define list of file names for Ra data
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

% Define time range for data extraction
years = 2015:2023; % Year range to extract
num_years = length(years);
num_months = 12;

% Dynamically read dimensions from first file
ra_sample = ncread(file_names{1}, 'ra');
[lon_size, lat_size, ~] = size(ra_sample);

% Initialize matrices for monthly and annual data storage
annual_ra_combined = nan(180, 360, num_months, num_years, length(file_names)); % Monthly data
annual_ra_sum = nan(180, 360, num_years, length(file_names));                 % Annual data

% Define days in each month (non-leap year as baseline)
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Process each file
for f = 1:length(file_names)
    f % Display current file number
    % Read 'ra' variable
    ra_data = ncread(file_names{f}, 'ra'); % Assume ra data is 3D: lon x lat x time
    [~, ~, time_dim] = size(ra_data);

    % Special handling for ISBA files
    if contains(file_names{f}, 'ISBA')
        % ISBA files: 360x150 -> 150x360, add NaN to pad to 180x360
        ra_data_rotated = rot90(ra_data, 1); % Rotate to 150x360
        ra_data_rotated = cat(1, ra_data_rotated, nan(30, 360, time_dim)); % Pad to 180x360
    elseif contains(file_names{f}, 'CLM6')
        % CLM6 files: special rotation and reorganization
        ra_data_rotated = rot90(ra_data, 1);
        ra_data_rotated = [ra_data_rotated(:,181:360,:), ra_data_rotated(:,1:180,:)];
    else
        % Regular files: 360x180 -> 180x360
        ra_data_rotated = rot90(ra_data, 1);
    end

    % Determine extraction time range
    start_idx = time_dim - num_years * num_months + 1; % Start from 2015
    end_idx = time_dim;                               % End at last month of 2023

    % Extract 2015-2023 data
    ra_sub_data = ra_data_rotated(:, :, start_idx:end_idx); % Extract last 9 years of data

    % Process data (monthly accumulation with leap year consideration)
    for y = 1:num_years
        current_year = years(y);

        % Check if leap year
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % Leap year: February has 29 days
        else
            days_in_month(2) = 28; % Non-leap year: February has 28 days
        end

        % Calculate seconds per month for current year
        seconds_per_month = days_in_month * 24 * 3600;

        % Extract current year's data
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        ra_yearly_data = ra_sub_data(:, :, month_start_idx:month_end_idx);

        % Convert units and accumulate
        ra_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            % Convert from kg/m²/s to gC/m²/month
            ra_monthly(:, :, m) = ra_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % Store in monthly matrix
        annual_ra_combined(:, :, :, y, f) = ra_monthly;

        % Sum 12 months of data and store in annual matrix
        annual_ra_sum(:, :, y, f) = sum(ra_monthly, 3); % Sum along 3rd dimension
    end
end

%% TRENDY Rh (Heterotrophic Respiration) Data Processing

% Define list of file names for Rh data
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

% Define time range for data extraction
years = 2015:2023; % Year range to extract
num_years = length(years);
num_months = 12;

% Dynamically read dimensions from first file
rh_sample = ncread(file_names{1}, 'rh');
[lon_size, lat_size, ~] = size(rh_sample);

% Initialize matrices for monthly and annual data storage
annual_rh_combined = nan(180, 360, num_months, num_years, length(file_names)); % Monthly data
annual_rh_sum = nan(180, 360, num_years, length(file_names));                 % Annual data

% Define days in each month (non-leap year as baseline)
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Process each file
for f = 1:length(file_names)
    f % Display current file number
    % Read 'rh' variable
    rh_data = ncread(file_names{f}, 'rh'); % Assume rh data is 3D: lon x lat x time
    [~, ~, time_dim] = size(rh_data);

    % Special handling for ISBA files
    if contains(file_names{f}, 'ISBA')
        % ISBA files: 360x150 -> 150x360, add NaN to pad to 180x360
        rh_data_rotated = rot90(rh_data, 1); % Rotate to 150x360
        rh_data_rotated = cat(1, rh_data_rotated, nan(30, 360, time_dim)); % Pad to 180x360
    elseif contains(file_names{f}, 'CLM6')
        % CLM6 files: special rotation and reorganization
        rh_data_rotated = rot90(rh_data, 1);
        rh_data_rotated = [rh_data_rotated(:,181:360,:), rh_data_rotated(:,1:180,:)];
    else
        % Regular files: 360x180 -> 180x360
        rh_data_rotated = rot90(rh_data, 1);
    end

    % Determine extraction time range
    start_idx = time_dim - num_years * num_months + 1; % Start from 2015
    end_idx = time_dim;                               % End at last month of 2023

    % Extract 2015-2023 data
    rh_sub_data = rh_data_rotated(:, :, start_idx:end_idx); % Extract last 9 years of data

    % Process data (monthly accumulation with leap year consideration)
    for y = 1:num_years
        current_year = years(y);

        % Check if leap year
        if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
            days_in_month(2) = 29; % Leap year: February has 29 days
        else
            days_in_month(2) = 28; % Non-leap year: February has 28 days
        end

        % Calculate seconds per month for current year
        seconds_per_month = days_in_month * 24 * 3600;

        % Extract current year's data
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        rh_yearly_data = rh_sub_data(:, :, month_start_idx:month_end_idx);

        % Convert units and accumulate
        rh_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            % Convert from kg/m²/s to gC/m²/month
            rh_monthly(:, :, m) = rh_yearly_data(:, :, m) * 1000 * seconds_per_month(m);
        end

        % Store in monthly matrix
        annual_rh_combined(:, :, :, y, f) = rh_monthly;

        % Sum 12 months of data and store in annual matrix
        annual_rh_sum(:, :, y, f) = sum(rh_monthly, 3); % Sum along 3rd dimension
    end
end

%% TRENDY tas (Temperature) Data Processing

% Define list of file names for temperature data
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

% Define time range for data extraction
years = 2015:2023; % Year range to extract
num_years = length(years);
num_months = 12;

% Dynamically read dimensions from first file
tas_sample = ncread(file_names{1}, 'tas');
[lon_size, lat_size, ~] = size(tas_sample);

% Initialize matrices for monthly and annual data storage
annual_tas_combined = nan(180, 360, num_months, num_years, length(file_names)); % Monthly data
annual_tas_sum = nan(180, 360, num_years, length(file_names));                 % Annual data

% Process each file
for f = 1:length(file_names)
    f % Display current file number
    % Read temperature variable (different variable name for ISBA)
    if contains(file_names{f}, 'ISBA')
        tas_data = ncread(file_names{f}, 'tasLut'); % Assume tas data is 3D: lon x lat x time
        [~, ~, time_dim] = size(tas_data);
    else
        tas_data = ncread(file_names{f}, 'tas'); % Assume tas data is 3D: lon x lat x time
        [~, ~, time_dim] = size(tas_data);
    end

    % Special handling for ISBA files
    if contains(file_names{f}, 'ISBA')
        % ISBA files: 360x150 -> 150x360, add NaN to pad to 180x360
        tas_data_rotated = rot90(tas_data, 1); % Rotate to 150x360
        tas_data_rotated = cat(1, tas_data_rotated, nan(30, 360, time_dim)); % Pad to 180x360
    elseif contains(file_names{f}, 'CLM6')
        % CLM6 files: special rotation and reorganization
        tas_data_rotated = rot90(tas_data, 1);
        tas_data_rotated = [tas_data_rotated(:,181:360,:), tas_data_rotated(:,1:180,:)];
    else
        % Regular files: 360x180 -> 180x360
        tas_data_rotated = rot90(tas_data, 1);
    end

    % Determine extraction time range
    start_idx = time_dim - num_years * num_months + 1; % Start from 2015
    end_idx = time_dim;                               % End at last month of 2023

    % Extract 2015-2023 data
    tas_sub_data = tas_data_rotated(:, :, start_idx:end_idx); % Extract last 9 years of data

    % Process data (convert from Kelvin to Celsius)
    for y = 1:num_years
        % Extract current year's data
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        tas_yearly_data = tas_sub_data(:, :, month_start_idx:month_end_idx);

        % Convert units (Kelvin to Celsius)
        tas_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            tas_monthly(:, :, m) = tas_yearly_data(:, :, m) - 273.15; % Convert K to °C
        end

        % Store in monthly matrix
        annual_tas_combined(:, :, :, y, f) = tas_monthly;

        % Calculate annual average temperature
        annual_tas_sum(:, :, y, f) = mean(tas_monthly, 3); % Average along 3rd dimension
    end
end

%% TRENDY mrso (Soil Moisture) Data Processing

% Define list of file names for soil moisture data
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

% Define time range for data extraction
years = 2015:2023; % Year range to extract
num_years = length(years);
num_months = 12;

% Dynamically read dimensions from first file
mrso_sample = ncread(file_names{1}, 'mrso');
[lon_size, lat_size, ~] = size(mrso_sample);

% Initialize matrices for monthly and annual data storage
annual_mrso_combined = nan(180, 360, num_months, num_years, length(file_names)); % Monthly data
annual_mrso_sum = nan(180, 360, num_years, length(file_names));                 % Annual data

% Process each file
for f = 1:length(file_names)
    f % Display current file number
    % Read 'mrso' variable
    mrso_data = ncread(file_names{f}, 'mrso'); % Assume mrso data is 3D: lon x lat x time
    [~, ~, time_dim] = size(mrso_data);

    % Special handling for ISBA files
    if contains(file_names{f}, 'ISBA')
        % ISBA files: 360x150 -> 150x360, add NaN to pad to 180x360
        mrso_data_rotated = rot90(mrso_data, 1); % Rotate to 150x360
        mrso_data_rotated = cat(1, mrso_data_rotated, nan(30, 360, time_dim)); % Pad to 180x360
    elseif contains(file_names{f}, 'CLM6')
        % CLM6 files: special rotation and reorganization
        mrso_data_rotated = rot90(mrso_data, 1);
        mrso_data_rotated = [mrso_data_rotated(:,181:360,:), mrso_data_rotated(:,1:180,:)];
    else
        % Regular files: 360x180 -> 180x360
        mrso_data_rotated = rot90(mrso_data, 1);
    end

    % Determine extraction time range
    start_idx = time_dim - num_years * num_months + 1; % Start from 2015
    end_idx = time_dim;                               % End at last month of 2023

    % Extract 2015-2023 data
    mrso_sub_data = mrso_data_rotated(:, :, start_idx:end_idx); % Extract last 9 years of data

    % Process soil moisture data (no unit conversion needed)
    for y = 1:num_years
        % Extract current year's data
        month_start_idx = (y - 1) * num_months + 1;
        month_end_idx = y * num_months;
        mrso_yearly_data = mrso_sub_data(:, :, month_start_idx:month_end_idx);

        % Store monthly data directly (no unit conversion)
        mrso_monthly = nan(180, 360, num_months);
        for m = 1:num_months
            mrso_monthly(:, :, m) = mrso_yearly_data(:, :, m);
        end

        % Store in monthly matrix
        annual_mrso_combined(:, :, :, y, f) = mrso_monthly;

        % Calculate annual average soil moisture
        annual_mrso_sum(:, :, y, f) = mean(mrso_monthly, 3); % Average along 3rd dimension
    end
end

% Save all processed variables to MATLAB file
save TRENDY_variable.mat annual_GPP_combined annual_GPP_sum annual_mrso_combined annual_mrso_sum annual_ra_combined annual_ra_sum annual_rh_combined ...
    annual_rh_sum annual_tas_combined annual_tas_sum
