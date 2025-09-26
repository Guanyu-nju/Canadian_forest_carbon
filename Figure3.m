%% Mask create
clc
clear
% Canada mask
pixel_mask=importdata('E:\phd_file\Boreal_North_America\Canada_mask.tif');
forest_mask=importdata("E:\phd_file\Boreal_North_America\region_lu.tif");
forest_mask(forest_mask~=1)=nan;
pixel_mask=pixel_mask.*forest_mask;

% Create environment
data = geotiffread('E:\phd_file\Boreal_North_America\region_lu.tif');
info = geotiffinfo('E:\phd_file\Boreal_North_America\region_lu.tif');
[m, n] = size(data);

% Extract latitude coordinates from the first column
k = 1;
for i = 1:m
    for j = 1:1
        [lat, lon] = pix2latlon(info.RefMatrix, i, j);   % Read latitude for all rows in first column
        lat_(k, :) = lat; % Store latitude data as a column
        k = k + 1;
    end
end

% Extract longitude coordinates from the first row
k = 1;
for ii = 1:1
    for jj = 1:n
        [lat, lon] = pix2latlon(info.RefMatrix, ii, jj);   % Read longitude for all columns in first row
        lon_(k, :) = lon;  % Store longitude data as a column
        k = k + 1;
    end
end

% Create coordinate grids using meshgrid
[lon1, lat1] = meshgrid(lon_, lat_);

Boundry = shaperead("E:\phd_file\yuling_shiliang\countries.shp");
bou_canX = [Boundry(:).X];
bou_canY= [Boundry(:).Y];

clearvars -except lon1 lat1 pixel_mask bou_canX bou_canY m n
%% 2023 anomaly

Temperature_year2023=importdata("E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\summer\Temperature_summer_2023.tif").*pixel_mask;
RZSM_year2023=importdata("E:\phd_file\Boreal_North_America\RZSM\GLDAS\summer\RZSM_summer_2023.tif").*pixel_mask;


for year=2015:2022

    TEM_year_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\summer\Temperature_summer_' num2str(year) '.tif']);
    TEM_year_mean(:,:,year-2014)=TEM_year_temp.*pixel_mask;

    RZSM_year_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\GLDAS\summer\RZSM_summer_' num2str(year) '.tif']);
    RZSM_year_mean(:,:,year-2014)=RZSM_year_temp.*pixel_mask;

end


TEM_year_mean=nanmean(TEM_year_mean,3);
TEM_absAnomaly_year=Temperature_year2023-TEM_year_mean;
RZSM_year_mean=nanmean(RZSM_year_mean,3);
RZSM_absAnomaly_year=RZSM_year2023-RZSM_year_mean;
%% 2023 summer zscore
%
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;


% time series create
for year=2015:2023

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\summer\ER_summer_' num2str(year) '.tif']);
    GCAS_ER_summer_list(year-2014)=nansum(nansum(ER_mean_temp.*area_grid))/(nansum(nansum(area_grid)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\summer\ER_summer_' num2str(year) '.tif']);
    GCB2024_ER_summer_list(year-2014)=nansum(nansum(ER_mean_temp.*area_grid))/(nansum(nansum(area_grid)));

    TEM_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\summer\Temperature_summer_' num2str(year) '.tif']);
    TEM_summer_list(year-2014)=nansum(nansum(TEM_temp.*area_grid))/(nansum(nansum(area_grid)));

    RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\GLDAS\summer\RZSM_summer_' num2str(year) '.tif']);
    RZSM_summer_list(year-2014)=nansum(nansum(RZSM_temp.*area_grid))/(nansum(nansum(area_grid)));

    GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_' num2str(year) '.tif']);
    GPP_summer_list(year-2014)=nansum(nansum(GPP_temp.*area_grid))/(nansum(nansum(area_grid)));

end

%%

ff=figure
% t = tiledlayout(4,6);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
set(gcf,'unit','pixels','position',[729,56,1246,1042]);

% part 1 *********************************************************************************
% n1=nexttile([2,3])
a=axes('Position',[0.050205457463884,0.661363475881282,0.426607812295437,0.361153986239765]);
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% fill shape
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,TEM_absAnomaly_year,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(a,nclCM(79));
caxis([-2,2]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-2:1:2]);
set(get(h,'ylabel'),'string','°C','fontsize',14);
title('Temperature anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.0895067224182208 1.0892839097071 0],'FontName','Arial','FontSize',20,'fontweight','bold')



% part 2 *********************************************************************************
b=axes('Position',[0.528378812500409,0.661363475881282,0.426607812295437,0.361153986239765]);

m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,RZSM_absAnomaly_year,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); 
colormap(b,flipud(nclCM(79)));
caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','kg m^{-2}','fontsize',14);
title('RZSM anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.0895067224182208 1.0892839097071 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 3 *********************************************************************************
% Create scatter plot for Temperature vs GCAS TER with quadratic fit
c = axes('Position', [0.064205457463884 0.40595009596929 0.255216693418941 0.26807634918777]);
scatter(TEM_summer_list(:), GCAS_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Quadratic polynomial fitting (2nd degree)
[p2, s] = polyfit(TEM_summer_list(:), GCAS_ER_summer_list(:), 2);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(TEM_summer_list(:)), max(TEM_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p2, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.0989794543340237 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14, 'Color', 'k');

% Highlight 2023 data point
text(TEM_summer_list(end), GCAS_ER_summer_list(end), '2023\rightarrow ', 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Calculate and display p-value using F-test
P_value = F_test(TEM_summer_list(:), GCAS_ER_summer_list(:), 2);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.0989794543340237 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Format the plot
box on
set(gca, 'XTick', [13:1:16], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [13,16]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
xlabel('Temperature (°C)', 'FontName', 'Arial', 'FontSize', 14);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
text('string', 'c', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');


% part 4 *********************************************************************************
% Create scatter plot for RZSM vs GCAS TER with quadratic fit
d = axes('Position', [0.392295345305248 0.40595009596929 0.255216693418941 0.26807634918777]);
scatter(RZSM_summer_list(:), GCAS_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Quadratic polynomial fitting (2nd degree)
[p1, s] = polyfit(RZSM_summer_list(:), GCAS_ER_summer_list(:), 2);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(RZSM_summer_list(:)), max(RZSM_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p1, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.63127519228908 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14);

% Calculate and display p-value using F-test
P_value = F_test(RZSM_summer_list(:), GCAS_ER_summer_list(:), 2);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.63127519228908 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Highlight 2023 data point
text(RZSM_summer_list(end), GCAS_ER_summer_list(end), ' \leftarrow2023', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Format the plot
box on
set(gca, 'XTick', [320:25:420], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [320,420]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
xlabel('RZSM (kg m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set x-axis label
text('string', 'd', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');


% part 5 *********************************************************************************
% Create scatter plot for GPP vs GCAS TER with linear fit
e = axes('Position', [0.726805778893001 0.40595009596929 0.255216693418941 0.26807634918777]);
scatter(GPP_summer_list(:), GCAS_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Linear fitting (1st degree)
[p2, s] = polyfit(GPP_summer_list(:), GCAS_ER_summer_list(:), 1);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(GPP_summer_list(:)), max(GPP_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p2, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.63127519228908 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14);

% Highlight 2023 data point
text(GPP_summer_list(end), GCAS_ER_summer_list(end), '2023\rightarrow ', 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Calculate and display p-value using F-test
P_value = F_test(GPP_summer_list(:), GCAS_ER_summer_list(:), 1);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.63127519228908 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Format the plot
box on
set(gca, 'XTick', [410:5:435], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [410,435]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
xlabel('GPP (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
text('string', 'e', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');


% part 6 *********************************************************************************
% Create scatter plot for Temperature vs GCB2024 TER with quadratic fit
f = axes('Position', [0.064205457463884 0.0580897248100616 0.255216693418941 0.26807634918777]);
scatter(TEM_summer_list(:), GCB2024_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Quadratic polynomial fitting (2nd degree)
[p2, s] = polyfit(TEM_summer_list(:), GCB2024_ER_summer_list(:), 2);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(TEM_summer_list(:)), max(TEM_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p2, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.0989794543340237 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14);

% Calculate and display p-value using F-test
P_value = F_test(TEM_summer_list(:), GCB2024_ER_summer_list(:), 2);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.0989794543340237 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Highlight 2023 data point
text(TEM_summer_list(end), GCB2024_ER_summer_list(end), '2023\rightarrow ', 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Format the plot
box on
set(gca, 'XTick', [13:1:16], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [13,16]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
xlabel('Temperature (°C)', 'FontName', 'Arial', 'FontSize', 14);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
text('string', 'f', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');


% part 7 *********************************************************************************
% Create scatter plot for RZSM vs GCB2024 TER with quadratic fit
g = axes('Position', [0.392295345305248 0.0580897248100616 0.255216693418941 0.26807634918777]);
scatter(RZSM_summer_list(:), GCB2024_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Quadratic polynomial fitting (2nd degree)
[p1, s] = polyfit(RZSM_summer_list(:), GCB2024_ER_summer_list(:), 2);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(RZSM_summer_list(:)), max(RZSM_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p1, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.63127519228908 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14);

% Calculate and display p-value using F-test
P_value = F_test(RZSM_summer_list(:), GCB2024_ER_summer_list(:), 2);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.63127519228908 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Highlight 2023 data point
text(RZSM_summer_list(end), GCB2024_ER_summer_list(end), ' \leftarrow2023', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Format the plot
box on
set(gca, 'XTick', [320:25:420], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [320,420]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
xlabel('RZSM (kg m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set x-axis label
text('string', 'g', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');


% part 8 *********************************************************************************
% Create scatter plot for GPP vs GCB2024 TER with linear fit
h = axes('Position', [0.726805778893001 0.0580897248100616 0.255216693418941 0.26807634918777]);
scatter(GPP_summer_list(:), GCB2024_ER_summer_list(:), 'filled', 'MarkerFaceColor', [173,87,38]/255, 'SizeData', 80); 
hold on

% Linear fitting (1st degree)
[p2, s] = polyfit(GPP_summer_list(:), GCB2024_ER_summer_list(:), 1);  % p contains fitting coefficients, s contains fitting statistics
x1 = linspace(min(GPP_summer_list(:)), max(GPP_summer_list(:)));
% Calculate fitted values
[y1, delta] = polyval(p2, x1, s);  % y1 contains fitted values, delta contains confidence interval deviation
plot(x1, y1, '--', 'LineWidth', 1.5, 'color', 'k'); 
hold on
% Plot confidence interval shading
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display R-squared value
R2 = s.rsquared;
R2 = roundn(R2, -2);
sig_text = ['\itR\rm^2 = ' num2str(R2)];
text('string', sig_text, 'Units', 'normalized', 'position', [0.63127519228908 0.166989437289201 0], 'FontName', 'Arial', 'FontSize', 14);

% Calculate and display p-value using F-test
P_value = F_test(GPP_summer_list(:), GCB2024_ER_summer_list(:), 1);
P_value = fix(P_value * 100) / 100;
P_text = ['\itP\rm = ' num2str(P_value)];
text('string', P_text, 'Units', 'normalized', 'position', [0.63127519228908 0.0687438232541132 0], 'FontName', 'Arial', 'FontSize', 14);

% Highlight 2023 data point
text(GPP_summer_list(end), GCB2024_ER_summer_list(end), '2023\rightarrow ', 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial', 'fontweight', 'bold', 'Color', 'k');

% Format the plot
box on
set(gca, 'XTick', [410:5:435], 'FontSize', 14, 'FontName', 'Arial', 'XLim', [410,435]);
set(gca, 'YTick', [230:40:350], 'FontSize', 14, 'FontName', 'Arial', 'YLim', [230,350]);
xlabel('GPP (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);
ylabel('TER (gC m^{-2})', 'FontName', 'Arial', 'FontSize', 14);  % Set y-axis label
text('string', 'h', 'Units', 'normalized', 'position', [-0.233835411577201 1.02242670669966 0], 'FontName', 'Arial', 'FontSize', 20, 'fontweight', 'bold');

% Save the figure
result = ['E:\phd_file\Boreal_North_America\Result\V9\anomaly2023_summer_scatter_map.png'];
% print(result, ff, '-r600', '-dpng');



