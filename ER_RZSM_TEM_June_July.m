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

area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;


% 时间序列均值
for year=2015:2023


    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_6.tif']);
    ER_June_list(year-2014)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_7.tif']);
    ER_July_list(year-2014)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_8.tif']);
    ER_Aug_list(year-2014)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    TEM_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_6.tif']);
    TEM_June_list(year-2014)=nansum(nansum(TEM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    TEM_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_7.tif']);
    TEM_July_list(year-2014)=nansum(nansum(TEM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    TEM_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_8.tif']);
    TEM_Aug_list(year-2014)=nansum(nansum(TEM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    ERA5_RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_6.tif']);
    ERA5_RZSM_June_list(year-2014)=nansum(nansum(ERA5_RZSM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    ERA5_RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_7.tif']);
    ERA5_RZSM_July_list(year-2014)=nansum(nansum(ERA5_RZSM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    ERA5_RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_8.tif']);
    ERA5_RZSM_Aug_list(year-2014)=nansum(nansum(ERA5_RZSM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

end

%%
f=figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% part 1 *********************************************************************************
nexttile
scatter(ERA5_RZSM_June_list(:),ER_June_list(:),[],TEM_June_list(:),'filled','SizeData',80); hold on

% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(ERA5_RZSM_June_list(:),ER_June_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(ERA5_RZSM_June_list(:)),max(ERA5_RZSM_June_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.62140648978107 0.0718794820929804 0],'FontName','Arial','FontSize',14,'Color','k')
text(ERA5_RZSM_June_list(end),ER_June_list(end),' \leftarrow2023','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)

colormap(flipud(hot));
caxis([11,15]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[11:2:15]);
set(get(h,'ylabel'),'string','June TEM (°C)','fontsize',14);
box on
ylim([0.23,0.5])
box on
ylabel('June ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('June RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[0.3:0.02:0.36],'FontSize',14,'FontName','Arial','XLim',[0.308,0.365]);
set(gca, 'YTick',[0.3:0.1:0.5],'FontSize',14,'FontName','Arial','YLim',[0.27,0.55]);
text('string','a','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 1 *********************************************************************************
nexttile
scatter(TEM_June_list(:),ER_June_list(:),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
scatter(TEM_June_list(end),ER_June_list(end),'filled','MarkerFaceColor',[0,107,172]/255,'SizeData',80); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(TEM_June_list(:),ER_June_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(TEM_June_list(:)),max(TEM_June_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.122876230565749 0.0680739307457945 0],'FontName','Arial','FontSize',14,'Color','k')
text(TEM_June_list(end),ER_June_list(end),'2023\rightarrow ','HorizontalAlignment','right','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)
box on
ylim([0.23,0.5])
box on
ylabel('June ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('June TEM (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[11:2:15],'FontSize',14,'FontName','Arial','XLim',[11,15]);
set(gca, 'YTick',[0.3:0.1:0.5],'FontSize',14,'FontName','Arial','YLim',[0.27,0.55]);
text('string','b','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')


% part 5 *********************************************************************************
nexttile
scatter(ERA5_RZSM_July_list(:),ER_July_list(:),[],TEM_July_list(:),'filled','SizeData',80); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(ERA5_RZSM_July_list(:),ER_July_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(ERA5_RZSM_July_list(:)),max(ERA5_RZSM_July_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
% R2 = sum((polyval(p,RZSM_summer_list(:)) - mean(ER_summer_list(:))).^2) / sum((ER_summer_list(:) - mean(ER_summer_list(:))).^2);  % 确定系数
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.62140648978107 0.0718794820929804 0],'FontName','Arial','FontSize',14,'Color','k')
text(ERA5_RZSM_July_list(end),ER_July_list(end),' \leftarrow2023','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)

colormap(flipud(hot));
caxis([15,17]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[15:1:17]);
set(get(h,'ylabel'),'string','July TEM (°C)','fontsize',14);
box on
ylabel('July ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('July RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[0.3:0.01:0.35],'FontSize',14,'FontName','Arial','XLim',[0.295,0.33]);
set(gca, 'YTick',[0.5:0.1:0.7],'FontSize',14,'FontName','Arial','YLim',[0.49,0.74]);
text('string','c','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 2 *********************************************************************************
nexttile
scatter(TEM_July_list(:),ER_July_list(:),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
scatter(TEM_July_list(end),ER_July_list(end),'filled','MarkerFaceColor',[0,107,172]/255,'SizeData',80); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(TEM_July_list(:),ER_July_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(TEM_July_list(:)),max(TEM_July_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.122876230565749 0.0680739307457945 0],'FontName','Arial','FontSize',14,'Color','k')
text(TEM_July_list(end),ER_July_list(end),'2023\rightarrow ','HorizontalAlignment','right','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)
box on
ylabel('July ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('July TEM (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[15:1:18],'FontSize',14,'FontName','Arial','XLim',[15,17.5]);
set(gca, 'YTick',[0.5:0.1:0.7],'FontSize',14,'FontName','Arial','YLim',[0.49,0.74]);
text('string','d','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 5 *********************************************************************************
nexttile
scatter(ERA5_RZSM_Aug_list(:),ER_Aug_list(:),[],TEM_Aug_list(:),'filled','SizeData',80); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(ERA5_RZSM_Aug_list(:),ER_Aug_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(ERA5_RZSM_Aug_list(:)),max(ERA5_RZSM_Aug_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
% R2 = sum((polyval(p,RZSM_summer_list(:)) - mean(ER_summer_list(:))).^2) / sum((ER_summer_list(:) - mean(ER_summer_list(:))).^2);  % 确定系数
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.62140648978107 0.0718794820929804 0],'FontName','Arial','FontSize',14,'Color','k')
text(ERA5_RZSM_Aug_list(end),ER_Aug_list(end),' \leftarrow2023','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)

colormap(flipud(hot));
caxis([14,16]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[14:1:16]);
set(get(h,'ylabel'),'string','August TEM (°C)','fontsize',14);
box on
ylabel('August ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('August RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[0.29:0.01:0.32],'FontSize',14,'FontName','Arial','XLim',[0.29,0.323]);
set(gca, 'YTick',[0.45:0.1:0.65],'FontSize',14,'FontName','Arial','YLim',[0.45,0.67]);
text('string','e','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 2 *********************************************************************************
nexttile
scatter(TEM_Aug_list(:),ER_Aug_list(:),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
scatter(TEM_Aug_list(end),ER_Aug_list(end),'filled','MarkerFaceColor',[0,107,172]/255,'SizeData',80); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(TEM_Aug_list(:),ER_Aug_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(TEM_Aug_list(:)),max(TEM_Aug_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.122876230565749 0.0680739307457945 0],'FontName','Arial','FontSize',14,'Color','k')
text(TEM_Aug_list(end),ER_Aug_list(end),'2023\rightarrow ','HorizontalAlignment','right','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)
box on
ylabel('August ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('August TEM (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[13:1:16],'FontSize',14,'FontName','Arial','XLim',[13,16.2]);
set(gca, 'YTick',[0.45:0.1:0.65],'FontSize',14,'FontName','Arial','YLim',[0.45,0.67]);
text('string','f','Units','normalized','position',[-0.187893000203482 1.02444976560022 0],'FontName','Arial','FontSize',14,'fontweight','bold')
set(gcf,'unit','centimeters','position',[12.250208333333338,7.366,22.484291666666664,25.971500000000006]);
result=['E:\phd_file\Boreal_North_America\Result\V6\ER_RZSM_TEM_June_July.png']
% print(result,f,'-r600','-dpng');
 
