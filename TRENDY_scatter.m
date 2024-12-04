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
%%
load TRENDY_variable.mat
annual_GPP_combined=annual_GPP_combined(:,:,:,:,1:end-1);
annual_mrso_combined=annual_mrso_combined(:,:,:,:,1:end-1);
annual_ra_combined=annual_ra_combined(:,:,:,:,1:end-1);
annual_rh_combined=annual_rh_combined(:,:,:,:,1:end-1);
annual_tas_combined=annual_tas_combined(:,:,:,:,1:end-1);
annual_GPP_sum=annual_GPP_sum(:,:,:,1:end-1);
annual_mrso_sum=annual_mrso_sum(:,:,:,1:end-1);
annual_ra_sum=annual_ra_sum(:,:,:,1:end-1);
annual_rh_sum=annual_rh_sum(:,:,:,1:end-1);
annual_tas_sum=annual_tas_sum(:,:,:,1:end-1);
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;

% 遍历每个文件计算ER NEP
for f = 1:size(annual_GPP_combined,5)
    for y = 1:size(annual_GPP_combined,4)
        % 计算 ER (逐月)
        annual_ER_combined(:, :, :, y, f) = annual_ra_combined(:, :, :, y, f) + annual_rh_combined(:, :, :, y, f);
    end
end

% 遍历每个文件
for f = 1:size(annual_ER_combined,5)
    for y = 1:size(annual_ER_combined,4)
        % 提取当前年的夏季数据 (6-8月)
        ER_summer = annual_ER_combined(:, :, 6:8, y, f); % 提取第 6-8 月的数据
        ER_summer=sum(ER_summer,3);
        mrso_summer = annual_mrso_combined(:, :, 6:8, y, f); % 提取第 6-8 月的数据
        mrso_summer=nanmean(mrso_summer,3);
        tas_summer = annual_tas_combined(:, :, 6:8, y, f); % 提取第 6-8 月的数据
        tas_summer=nanmean(tas_summer,3);

        % 按面积加权求和，得到单个文件的夏季总量
        summer_ER_list(f, y) = nansum(nansum(ER_summer .* area_grid))/(10^15);
        summer_mrso_list(f, y) = nansum(nansum(mrso_summer.*area_grid))/(nansum(nansum(area_grid)));
        summer_tas_list(f, y) = nansum(nansum(tas_summer.*area_grid))/(nansum(nansum(area_grid)));
    end
end
%%
load CMIP6_summer_variable_result.mat
ff=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
s1=scatter(SM_list2,ER_list2,'filled','MarkerFaceColor',[68,133,199]/255,'SizeData',80); hold on; box on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(SM_list2,ER_list2,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(SM_list2(:)),max(SM_list2(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.0752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
set(gca,'FontSize',14,'FontName','Arial');
xlim([1348,1374])
set(gca, 'XTick',[1350:10:1370],'FontSize',14,'FontName','Arial');
ylim([2.2,2.45])
xlabel('RZSM (kg m^{-2})')
ylabel('ER (PgC)')
text('string','a','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',14,'fontweight','bold')

nexttile
scatter(TEM_list2,ER_list2,'filled','MarkerFaceColor',[68,133,199]/255,'SizeData',80); hold on; box on
% title(filelist{i})
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(TEM_list2,ER_list2,2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(TEM_list2(:)),max(TEM_list2(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.0752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
set(gca,'FontSize',14,'FontName','Arial');
xlim([13.7,15.2])
ylim([2.2,2.45])
xlabel('TEM (°C)')
ylabel('ER (PgC)')
text('string','b','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',14,'fontweight','bold')

nexttile
s2=scatter(mean(summer_mrso_list),mean(summer_ER_list),'filled','MarkerFaceColor',[132,186,66]/255,'SizeData',80); hold on; box on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(mean(summer_mrso_list),mean(summer_ER_list),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(mean(summer_mrso_list)),max(mean(summer_mrso_list)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.0752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
xlim([1005,1027])
set(gca, 'XTick',[1007:6:1025],'FontSize',14,'FontName','Arial');

xlabel('RZSM (kg m^{-2})')
ylabel('ER (PgC)')
text('string','c','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial');

nexttile
scatter(mean(summer_tas_list),mean(summer_ER_list),'filled','MarkerFaceColor',[132,186,66]/255,'SizeData',80); hold on; box on
[p,s] = polyfit(mean(summer_tas_list),mean(summer_ER_list),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(mean(summer_tas_list)),max(mean(summer_tas_list)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.68839347194506 0.0752941473522928 0],'FontName','Arial','FontSize',14,'Color','k')
xlim([12.8,14.2])

xlabel('TEM (°C)','FontName','Arial','FontSize',14)
ylabel('ER (PgC)','FontName','Arial','FontSize',14)
text('string','d','Units','normalized','position',[-0.190626006158388 1.04594302254788 0],'FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial');
lgd = legend([s1 s2],{'CMIP6','TRENDY'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','on');
lgd.Layout.Tile = 'east';
lgd.ItemTokenSize = [15 15];
set(gcf,'unit','centimeters','position',[26.431875000000005,13.176250000000001,24.421041666666667,19.552708333333335]);
result=['E:\phd_file\Boreal_North_America\Result\V6\TRENDY_CMIP6_scatter.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCN","ORCHIDEE","SDGVM","VISIT","VISIT-UT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ii=1;
for i=1:8
    nexttile
    scatter(summer_mrso_list(i,:),summer_ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    % 使用二项式拟合 (二次多项式)
    [p,s] = polyfit(summer_mrso_list(i,:),summer_ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(summer_mrso_list(i,:)),max(summer_mrso_list(i,:)));
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
    scatter(summer_tas_list(i,:),summer_ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(summer_tas_list(i,:),summer_ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(summer_tas_list(i,:)),max(summer_tas_list(i,:)));
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
padding_factor = 0.25;  % Adjust padding as needed

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
result=['E:\phd_file\Boreal_North_America\Result\V6\TRENDY_indiviual_scatter1.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCN","ORCHIDEE","SDGVM","VISIT","VISIT-UT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x"};
ff=figure
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ii=1;
for i=9:14
    nexttile
    scatter(summer_mrso_list(i,:),summer_ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    % 使用二项式拟合 (二次多项式)
    [p,s] = polyfit(summer_mrso_list(i,:),summer_ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(summer_mrso_list(i,:)),max(summer_mrso_list(i,:)));
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
    scatter(summer_tas_list(i,:),summer_ER_list(i,:),'filled','MarkerFaceColor',[173,87,38]/255); hold on; box on
    [p,s] = polyfit(summer_tas_list(i,:),summer_ER_list(i,:),2);  % p 是拟合系数，S 是拟合统计信息
    x1=linspace(min(summer_tas_list(i,:)),max(summer_tas_list(i,:)));
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
padding_factor = 0.25;  % Adjust padding as needed

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
result=['E:\phd_file\Boreal_North_America\Result\V6\TRENDY_indiviual_scatter2.png']
% print(result,ff,'-r600','-dpng' );

