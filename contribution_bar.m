%% Mask creation
clc
clear
% Canada mask
pixel_mask = importdata('E:\phd_file\Boreal_North_America\Canada_mask.tif');
forest_mask = importdata("E:\phd_file\Boreal_North_America\region_lu.tif");
forest_mask(forest_mask ~= 1) = nan;
pixel_mask = pixel_mask .* forest_mask;

% Create coordinate environment
data = geotiffread('E:\phd_file\Boreal_North_America\region_lu.tif');
info = geotiffinfo('E:\phd_file\Boreal_North_America\region_lu.tif');
[m, n] = size(data);

% Extract latitude coordinates (first column)
k = 1;
for i = 1:m
    for j = 1:1
        [lat, lon] = pix2latlon(info.RefMatrix, i, j);   % Read latitude for all rows in first column
        lat_(k, :) = lat; % Store latitude data as a column
        k = k + 1;
    end
end

% Extract longitude coordinates (first row)
k = 1;
for ii = 1:1
    for jj = 1:n
        [lat, lon] = pix2latlon(info.RefMatrix, ii, jj);   % Read longitude for all columns in first row
        lon_(k, :) = lon;  % Store longitude data as a column
        k = k + 1;
    end
end

% Create coordinate grids
[lon1, lat1] = meshgrid(lon_, lat_);

% Load country boundary shapefile
Boundry = shaperead("E:\phd_file\yuling_shiliang\countries.shp");
bou_canX = [Boundry(:).X];
bou_canY = [Boundry(:).Y];

% Clean up workspace, keep only necessary variables
clearvars -except lon1 lat1 pixel_mask bou_canX bou_canY m n

%%
% Load area grid data and apply mask
area_grid = importdata("E:\phd_file\Boreal_North_America\degree2meter.tif") * 1000000 .* pixel_mask;

% Time series mean calculation
for year = 2015:2023
    % GCAS ER data
    ER_mean_temp = importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\summer\ER_summer_' num2str(year) '.tif']);
    GCAS_ER_list(year-2014, 1) = nansum(nansum(ER_mean_temp .* area_grid)) / (nansum(nansum(area_grid)));

    % GCB2024 ER data
    ER_mean_temp = importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\summer\ER_summer_' num2str(year) '.tif']);
    GCB_ER_list(year-2014, 1) = nansum(nansum(ER_mean_temp .* area_grid)) / (nansum(nansum(area_grid)));

    % Temperature data
    TEM_temp = importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\summer\Temperature_summer_' num2str(year) '.tif']);
    TEM_list(year-2014, 1) = nansum(nansum(TEM_temp .* area_grid)) / (nansum(nansum(area_grid)));

    % Root zone soil moisture data
    RZSM_temp = importdata(['E:\phd_file\Boreal_North_America\RZSM\GLDAS\summer\RZSM_summer_' num2str(year) '.tif']);
    RZSM_list(year-2014, 1) = nansum(nansum(RZSM_temp .* area_grid)) / (nansum(nansum(area_grid)));

    % GPP data
    GPP_temp = importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_' num2str(year) '.tif']);
    GPP_list(year-2014, 1) = nansum(nansum(GPP_temp .* area_grid)) / (nansum(nansum(area_grid)));
end

% Normalization (Z-score)
RZSM_norm = normalize(RZSM_list, 'range');
TEM_norm = normalize(TEM_list, 'range');
GPP_norm = normalize(GPP_list, 'range');

% Prepare feature matrices with polynomial terms
% Format: [RZSM^2, TEMP^2, RZSM, TEMP, GPP, intercept]
GCAS_X = [RZSM_norm.^2, TEM_norm.^2, RZSM_norm, TEM_norm, GPP_norm, ones(size(GCAS_ER_list))];
GCB2024_X = [RZSM_norm.^2, TEM_norm.^2, RZSM_norm, TEM_norm, GPP_norm, ones(size(GCB_ER_list))];

% Fit linear regression model for GCAS data
GCAS_beta = GCAS_X \ GCAS_ER_list;
% Calculate predicted values
GCAS_y_pred = GCAS_X * GCAS_beta;

% Fit linear regression model for GCB2024 data
GCB2024_beta = GCB2024_X \ GCB_ER_list;
% Calculate predicted values
GCB2024_y_pred = GCB2024_X * GCB2024_beta;

% % Calculate R² and RMSE (commented out)
% SS_total = sum((y - mean(y)).^2);
% SS_residual = sum((y - GCAS_y_pred).^2);
% R2 = 1 - (SS_residual / SS_total);

GCAS_Climate_beta=GCAS_beta;GCAS_Climate_beta(5)=0
GCAS_ER_climate_predict_result= GCAS_X * GCAS_Climate_beta;
GCAS_ER_predict_2023anomalies=GCAS_y_pred(end)-nanmean(GCAS_y_pred(1:end-1))
GCAS_ER_climate_predict_2023anomalies=GCAS_ER_climate_predict_result(end)-nanmean(GCAS_ER_climate_predict_result(1:end-1))
GCAS_result=[GCAS_ER_climate_predict_2023anomalies,GCAS_ER_predict_2023anomalies-GCAS_ER_climate_predict_2023anomalies];

GCB2024_Climate_beta=GCB2024_beta;GCB2024_Climate_beta(5)=0
GCB2024_ER_climate_predict_result= GCB2024_X * GCB2024_Climate_beta;
GCB2024_ER_predict_2023anomalies=GCB2024_y_pred(end)-nanmean(GCB2024_y_pred(1:end-1))
GCB2024_ER_climate_predict_2023anomalies=GCB2024_ER_climate_predict_result(end)-nanmean(GCB2024_ER_climate_predict_result(1:end-1))
GCB2024_result=[GCB2024_ER_climate_predict_2023anomalies,GCB2024_ER_predict_2023anomalies-GCB2024_ER_climate_predict_2023anomalies];

Total_result=[GCAS_result;GCB2024_result];

%%
ff=figure
set(gcf,'unit','pixels','position',[1000,877,855,361]);
b=bar(Total_result)
b(1).FaceColor = [67 138 200]/255;
b(2).FaceColor = [254 204 81]/255;
legend({'RZSM + temperature','GPP'},'FontName','Arial','FontSize',14,'Location','eastoutside')
ylabel('TER anomalies (gC m^{-2})','FontName','Arial','FontSize',14);
set(gca,'XTickLabel',{'GCAS-extra','GCB2024-satellite'},'FontName','Arial','fontsize',14)
result=['E:\phd_file\Boreal_North_America\Result\V9\Contribution_bar.png']
% print(result,ff,'-r600','-dpng');

