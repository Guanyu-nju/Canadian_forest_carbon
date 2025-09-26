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
%%
% Load TRENDY model variables
load TRENDY_variable.mat

% Load and prepare area grid data for weighted calculations
% Convert degree to meter and scale, then apply pixel mask
area_grid = importdata("E:\phd_file\Boreal_North_America\degree2meter.tif") * 1000000 .* pixel_mask;

% Calculate ER (Ecosystem Respiration) and NEP for each file
% Loop through each model file
for f = 1:size(annual_GPP_combined, 5)
    % Loop through each year
    for y = 1:size(annual_GPP_combined, 4)
        % Calculate monthly ER (Ecosystem Respiration = autotrophic respiration + heterotrophic respiration)
        annual_ER_combined(:, :, :, y, f) = annual_ra_combined(:, :, :, y, f) + annual_rh_combined(:, :, :, y, f);
    end
end

% Process summer data for each model file and year
% Loop through each model file
for f = 1:size(annual_ER_combined, 5)
    % Loop through each year
    for y = 1:size(annual_ER_combined, 4)
        % Extract summer data (June-August, months 6-8) for current year
        % Ecosystem Respiration summer data
        ER_summer = annual_ER_combined(:, :, 6:8, y, f); % Extract June-August data
        ER_summer = sum(ER_summer, 3); % Sum over summer months
        
        % GPP summer data
        GPP_summer = annual_GPP_combined(:, :, 6:8, y, f); % Extract June-August data
        GPP_summer = sum(GPP_summer, 3); % Sum over summer months
        
        % Soil moisture summer data (mean instead of sum)
        mrso_summer = annual_mrso_combined(:, :, 6:8, y, f); % Extract June-August data
        mrso_summer = nanmean(mrso_summer, 3); % Average over summer months
        
        % Temperature summer data (mean instead of sum)
        tas_summer = annual_tas_combined(:, :, 6:8, y, f); % Extract June-August data
        tas_summer = nanmean(tas_summer, 3); % Average over summer months

        % Calculate area-weighted summer totals/averages for each variable
        % Ecosystem Respiration summer total (area-weighted)
        summer_ER_list(f, y) = nansum(nansum(ER_summer .* area_grid)) / (nansum(nansum(area_grid)));
        
        % GPP summer total (area-weighted)
        summer_GPP_list(f, y) = nansum(nansum(GPP_summer .* area_grid)) / (nansum(nansum(area_grid)));
        
        % Soil moisture summer average (area-weighted)
        summer_mrso_list(f, y) = nansum(nansum(mrso_summer .* area_grid)) / (nansum(nansum(area_grid)));
        
        % Temperature summer average (area-weighted)
        summer_tas_list(f, y) = nansum(nansum(tas_summer .* area_grid)) / (nansum(nansum(area_grid)));
    end
end
summer_NEP_list=summer_GPP_list-summer_ER_list;

for i=1:length(summer_ER_list)

    summer_ER_anomaly_list(i)=summer_ER_list(i,end)-nanmean(summer_ER_list(i,1:end-1));
    summer_GPP_anomaly_list(i)=summer_GPP_list(i,end)-nanmean(summer_GPP_list(i,1:end-1));
    summer_NEP_anomaly_list(i)=summer_NEP_list(i,end)-nanmean(summer_NEP_list(i,1:end-1));
    
end
%%
ff=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[1000,529,780,709]);

nexttile
result1=[summer_NEP_anomaly_list;summer_GPP_anomaly_list;summer_ER_anomaly_list]'
b=boxchart(result1,'BoxFaceColor',[68,133,199]/255); hold on
box on
b.MarkerStyle = 'none';
set(gca,'XTickLabel',{'NEP','GPP','TER'},'FontName','Arial','fontsize',14)
ylabel('Flux anomalies (gC m^{-2})','FontName','Arial','FontSize',14)
text('string','a','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',18,'fontweight','bold')
hh=yline(0, '--k', 'LineWidth', 1); 
ylim([-24,40])


nexttile
scatter(mean(summer_tas_list),mean(summer_ER_list),'filled','MarkerFaceColor',[68,133,199]/255,'SizeData',80); hold on; box on
text(mean(summer_tas_list(:,end)),mean(summer_ER_list(:,end)),'2023\rightarrow ','HorizontalAlignment','right','FontSize',14,'FontName','Arial','fontweight','bold','Color', 'k')
[p,s] = polyfit(mean(summer_tas_list),mean(summer_ER_list),2);  
x1=linspace(min(mean(summer_tas_list)),max(mean(summer_tas_list)));
[y1, delta] = polyval(p, x1, s); 
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.1752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
P_value=F_test(mean(summer_tas_list),mean(summer_ER_list),2);
P_value=fix(P_value * 100) / 100;
if P_value > 0.01
    P_text=['\itP\rm = ' num2str(P_value)];
else
    P_text=['\itP\rm < 0.01'];
end
text('string',P_text,'Units','normalized','position',[0.68839347194506 0.0650481682454793 0],'FontName','Arial','FontSize',14,'Color','k')
% xlim([12.8,14.2])
xlabel('Temperature (°C)','FontName','Arial','FontSize',14)
ylabel('TER (gC m^{-2})','FontName','Arial','FontSize',14)
text('string','b','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',18,'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial');
ylim([440,475])
set(gca, 'YTick',[440:10:480],'FontSize',14,'FontName','Arial');

nexttile
scatter(mean(summer_mrso_list),mean(summer_ER_list),'filled','MarkerFaceColor',[68,133,199]/255,'SizeData',80); hold on; box on
text(mean(summer_mrso_list(:,end)),mean(summer_ER_list(:,end)),' \leftarrow2023','HorizontalAlignment','left','FontSize',14,'FontName','Arial','fontweight','bold','Color', 'k')
[p,s] = polyfit(mean(summer_mrso_list),mean(summer_ER_list),2); 
x1=linspace(min(mean(summer_mrso_list)),max(mean(summer_mrso_list)));
[y1, delta] = polyval(p, x1, s); 
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.0737755317457244 0.175294147352293 0],'FontName','Arial','FontSize',14,'Color','k')
P_value=F_test(mean(summer_mrso_list),mean(summer_ER_list),2);
P_value=fix(P_value * 100) / 100;
if P_value > 0.01
    P_text=['\itP\rm = ' num2str(P_value)];
else
    P_text=['\itP\rm < 0.01'];
end
text('string',P_text,'Units','normalized','position',[0.0737755317457244 0.0650481682454793 0],'FontName','Arial','FontSize',14,'Color','k')
xlim([1040,1065])
set(gca, 'XTick',[1040:10:1060,1065],'FontSize',14,'FontName','Arial');
ylim([440,475])
set(gca, 'YTick',[440:10:480],'FontSize',14,'FontName','Arial');
xlabel('RZSM (kg m^{-2})')
ylabel('TER (gC m^{-2})')
text('string','c','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',18,'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial');


nexttile
scatter(mean(summer_GPP_list),mean(summer_ER_list),'filled','MarkerFaceColor',[68,133,199]/255,'SizeData',80); hold on; box on
text(mean(summer_GPP_list(:,end)),mean(summer_ER_list(:,end)),'2023\rightarrow ','HorizontalAlignment','right','FontSize',14,'FontName','Arial','fontweight','bold','Color', 'k')
[p,s] = polyfit(mean(summer_GPP_list),mean(summer_ER_list),1); 
x1=linspace(min(mean(summer_GPP_list)),max(mean(summer_GPP_list)));
[y1, delta] = polyval(p, x1, s);
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['\itR\rm^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.1752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
P_value=F_test(mean(summer_GPP_list),mean(summer_ER_list),2);
P_value=fix(P_value * 100) / 100;
if P_value > 0.01
    P_text=['\itP\rm = ' num2str(P_value)];
else
    P_text=['\itP\rm < 0.01'];
end
text('string',P_text,'Units','normalized','position',[0.68839347194506 0.0650481682454793 0],'FontName','Arial','FontSize',14,'Color','k')
xlim([574,606])
ylim([440,475])
set(gca, 'YTick',[440:10:480],'FontSize',14,'FontName','Arial');

xlabel('GPP (gC m^{-2})','FontName','Arial','FontSize',14)
ylabel('TER (gC m^{-2})','FontName','Arial','FontSize',14)
text('string','d','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',18,'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial');



result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_scatter.png']
% print(result,ff,'-r600','-dpng' );

%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ii=1;
for i=1:14
    nexttile
    scatter(summer_mrso_list(i,end),summer_ER_list(i,end),'filled','MarkerFaceColor',[55,140,230]/255);hold on;
    scatter(summer_mrso_list(i,1:end-1),summer_ER_list(i,1:end-1),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    %
    [p,s] = polyfit(summer_mrso_list(i,:),summer_ER_list(i,:),2);
    x1=linspace(min(summer_mrso_list(i,:)),max(summer_mrso_list(i,:)));
    %
    [y1, delta] = polyval(p, x1, s);
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    %
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['\itR\rm^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.172310255575923 0],'FontName','Arial','FontSize',12,'Color','k')
    P_value=F_test(summer_mrso_list(i,:),summer_ER_list(i,:),2);
    P_value=fix(P_value * 100) / 100;
    if P_value > 0.01
        P_text=['\itP\rm = ' num2str(P_value)];
    else
        P_text=['\itP\rm < 0.01'];
    end
    text('string',P_text,'Units','normalized','position',[0.597736439528456 0.0650481682454793 0],'FontName','Arial','FontSize',12)
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('RZSM (kg m^{-2})','FontName','Arial','FontSize',12)
    ylabel('TER (gC m^{-2})','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',16,'fontweight','bold')
    ii=ii+1;
    % text(summer_mrso_list(i,end),summer_ER_list(i,end),' \leftarrow2023','FontSize',12,'FontName','Arial','fontweight','bold','Color', 'k')
end

% Assuming your subplots have already been created in a figure
figure_handle = gcf;  % Get current figure handle
axes_handles = findall(figure_handle, 'type', 'axes');  % Get all axes handles

% Define a padding factor to add extra space around scatter points
padding_factor = 0.65;  % Adjust padding as needed

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

        xlim(ax, [min(xData) - padding_factor*1.5 * xRange, max(xData) + padding_factor * xRange]);
        ylim(ax, [min(yData) - padding_factor * yRange, max(yData) + padding_factor*1.5 * yRange]);
    end
end

set(gcf,'unit','centimeters','position',[12.832291666666668,5.979583333333334,30.718125,28.86604166666668]);
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_indiviual_scatter1.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ii=1;
for i=1:14

    nexttile
    scatter(summer_tas_list(i,end),summer_ER_list(i,end),'filled','MarkerFaceColor',[55,140,230]/255);hold on;
    scatter(summer_tas_list(i,1:end-1),summer_ER_list(i,1:end-1),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(summer_tas_list(i,:),summer_ER_list(i,:),2);
    x1=linspace(min(summer_tas_list(i,:)),max(summer_tas_list(i,:)));
    %
    [y1, delta] = polyval(p, x1, s);
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['\itR\rm^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.172310255575923 0],'FontName','Arial','FontSize',12,'Color','k')
    P_value=F_test(summer_tas_list(i,:),summer_ER_list(i,:),2);
    P_value=fix(P_value * 100) / 100;
    if P_value > 0.01
        P_text=['\itP\rm = ' num2str(P_value)];
    else
        P_text=['\itP\rm < 0.01'];
    end
    text('string',P_text,'Units','normalized','position',[0.597736439528456 0.0650481682454793 0],'FontName','Arial','FontSize',12)
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('Temperature (°C)','FontName','Arial','FontSize',12)
    ylabel('TER (gC m^{-2})','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',16,'fontweight','bold')
    ii=ii+1;
    % text(summer_tas_list(i,end),summer_ER_list(i,end),' \leftarrow2023','FontSize',12,'FontName','Arial','fontweight','bold','Color', 'k')
end

% Assuming your subplots have already been created in a figure
figure_handle = gcf;  % Get current figure handle
axes_handles = findall(figure_handle, 'type', 'axes');  % Get all axes handles

% Define a padding factor to add extra space around scatter points
padding_factor = 0.3;  % Adjust padding as needed

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

set(gcf,'unit','centimeters','position',[12.832291666666668,5.979583333333334,30.718125,28.86604166666668]);
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_indiviual_scatter2.png']
% print(result,ff,'-r600','-dpng' );

%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ii=1;
for i=1:14

    nexttile
    scatter(summer_GPP_list(i,end),summer_ER_list(i,end),'filled','MarkerFaceColor',[55,140,230]/255);hold on;
    scatter(summer_GPP_list(i,1:end-1),summer_ER_list(i,1:end-1),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(summer_GPP_list(i,:),summer_ER_list(i,:),1);
    x1=linspace(min(summer_GPP_list(i,:)),max(summer_GPP_list(i,:)));
   
    [y1, delta] = polyval(p, x1, s);
    plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on

    fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
        'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    R2=s.rsquared;
    R2=roundn(R2,-2);
    sig_text=['\itR\rm^2 = ' num2str(R2)];
    text('string',sig_text,'Units','normalized','position',[0.597736439528456 0.172310255575923 0],'FontName','Arial','FontSize',12,'Color','k')
    P_value=F_test(summer_GPP_list(i,:),summer_ER_list(i,:),1);
    P_value=fix(P_value * 100) / 100;
    if P_value > 0.01
        P_text=['\itP\rm = ' num2str(P_value)];
    else
        P_text=['\itP\rm < 0.01'];
    end
    text('string',P_text,'Units','normalized','position',[0.597736439528456 0.0650481682454793 0],'FontName','Arial','FontSize',12)
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')
    xlabel('GPP (gC m^{-2})','FontName','Arial','FontSize',12)
    ylabel('TER (gC m^{-2})','FontName','Arial','FontSize',12)
    set(gca,'FontSize',12,'FontName','Arial');
    text('string',panellist{ii},'Units','normalized','position',[-0.202009629002369 1.13342409859307 0],'FontName','Arial','FontSize',16,'fontweight','bold')
    ii=ii+1;
    % text(summer_tas_list(i,end),summer_ER_list(i,end),' \leftarrow2023','FontSize',12,'FontName','Arial','fontweight','bold','Color', 'k')
end

% Assuming your subplots have already been created in a figure
figure_handle = gcf;  % Get current figure handle
axes_handles = findall(figure_handle, 'type', 'axes');  % Get all axes handles

% Define a padding factor to add extra space around scatter points
padding_factor = 0.3;  % Adjust padding as needed

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

set(gcf,'unit','centimeters','position',[12.832291666666668,5.979583333333334,30.718125,28.86604166666668]);
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_indiviual_scatter3.png']
% print(result,ff,'-r600','-dpng' );


