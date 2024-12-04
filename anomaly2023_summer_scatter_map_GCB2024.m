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

%% 2023 summer zscore
%
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;


% 时间序列均值
for year=2015:2023

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\summer\ER_summer_' num2str(year) '.tif']);
    ER_summer_list(year-2014)=nansum(nansum(ER_mean_temp.*area_grid/(10^15)));

    TEM_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\summer\Temperature_summer_' num2str(year) '.tif']);
    TEM_summer_list(year-2014)=nansum(nansum(TEM_temp.*area_grid))/(nansum(nansum(area_grid)));

    RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\summer\RZSM_summer_' num2str(year) '.tif']);
    RZSM_summer_list(year-2014)=nansum(nansum(RZSM_temp.*area_grid))/(nansum(nansum(area_grid)));

end

%%

f=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';


% part 3 *********************************************************************************
n3=nexttile
scatter(RZSM_summer_list(:),ER_summer_list(:),[],TEM_summer_list(:),'filled','SizeData',80); hold on

% 使用二项式拟合 (二次多项式)
[p1,s] = polyfit(RZSM_summer_list(:),ER_summer_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(RZSM_summer_list(:)),max(RZSM_summer_list(:)));
% 计算拟合值
[y1, delta] = polyval(p1, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
% R2 = sum((polyval(p,RZSM_summer_list(:)) - mean(ER_summer_list(:))).^2) / sum((ER_summer_list(:) - mean(ER_summer_list(:))).^2);  % 确定系数
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.209289492090949 0.16605008343821 0],'FontName','Arial','FontSize',14)
text(RZSM_summer_list(end),ER_summer_list(end),' \leftarrow2023','FontSize',14,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)

colormap(n3,flipud(hot));
caxis([13,16]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[13:1:16]);
set(get(h,'ylabel'),'string','Summer TEM (°C)','fontsize',14);
box on
ylabel('Summer ER (PgC)','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Summer RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[0.31:0.01:0.35],'FontSize',14,'FontName','Arial','XLim',[0.303,0.337]);
set(gca, 'YTick',[1.3:0.3:1.9],'FontSize',14,'FontName','Arial','YLim',[1.15,1.9]);
text('string','a','Units','normalized','position',[-0.16719463246739 1.03106688869709 0],'FontName','Arial','FontSize',14,'fontweight','bold')

[p2,s] = polyfit(TEM_summer_list(:),ER_summer_list(:),2);  % p 是拟合系数，S 是拟合统计信息


% part 1 *********************************************************************************
nexttile
[p3,s] = polyfit(TEM_summer_list(:),RZSM_summer_list(:),1);  % p 是拟合系数，S 是拟合统计信息
scatter(RZSM_summer_list(:)-polyval(p3,TEM_summer_list(:)),ER_summer_list(:)-polyval(p2,TEM_summer_list(:)),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
[p,s] = polyfit(RZSM_summer_list(:)-polyval(p3,TEM_summer_list(:)),ER_summer_list(:)-polyval(p2,TEM_summer_list(:)),2);  %多项式曲线拟合
x1=linspace(min(RZSM_summer_list(:)-polyval(p3,TEM_summer_list(:))),max(RZSM_summer_list(:)-polyval(p3,TEM_summer_list(:))));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.757206008913333 0.0690245876615243 0],'FontName','Arial','FontSize',14)
box on
ylabel('Summer ER|Summer TEM','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Summer RZSM|Summer TEM','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[-0.01:0.005:0.015],'FontSize',14,'FontName','Arial','XLim',[-0.01,0.016]);
set(gca, 'YTick',[-0.16:0.08:0.24],'FontSize',14,'FontName','Arial','YLim',[-0.18,0.18]);
text('string','b','Units','normalized','position',[-0.16719463246739 1.03106688869709 0],'FontName','Arial','FontSize',14,'fontweight','bold')


% part 4 *********************************************************************************
nexttile
[p4,s] = polyfit(RZSM_summer_list(:),TEM_summer_list(:),1);  % p 是拟合系数，S 是拟合统计信息
scatter(TEM_summer_list(:)-polyval(p4,RZSM_summer_list(:)),ER_summer_list(:)-polyval(p1,RZSM_summer_list(:)),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',80); hold on
[p,s] = polyfit(TEM_summer_list(:)-polyval(p4,RZSM_summer_list(:)),ER_summer_list(:)-polyval(p1,RZSM_summer_list(:)),2);  %多项式曲线拟合
x1=linspace(min(TEM_summer_list(:)-polyval(p4,RZSM_summer_list(:))),max(TEM_summer_list(:)-polyval(p4,RZSM_summer_list(:))));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.757206008913333 0.0690245876615243 0],'FontName','Arial','FontSize',14)
box on
ylabel('Summer ER|Summer RZSM','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Summer TEM|Summer RZSM','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[-1.2:0.4:0.8],'FontSize',14,'FontName','Arial','XLim',[-1.1,0.8]);
set(gca, 'YTick',[-0.18:0.06:0.12],'FontSize',14,'FontName','Arial','YLim',[-0.195,0.14]);
text('string','c','Units','normalized','position',[-0.16719463246739 1.03106688869709 0],'FontName','Arial','FontSize',14,'fontweight','bold')


n4=nexttile
[xData, yData, zData] = prepareSurfaceData( RZSM_summer_list, TEM_summer_list, ER_summer_list );
% 设置 fittype 和选项。
ft = fittype( 'lowess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 1;

% 对数据进行模型拟合。
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
% 绘制数据拟合图。
plot(fitresult); hold on
R2=gof.rsquare;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.757206008913333 0.364222311293267 0],'FontName','Arial','FontSize',14)
set(gca, 'ZTick',[1.3:0.3:1.9],'FontSize',14,'FontName','Arial','ZLim',[1.3,1.9]);
set(gca, 'YTick',[13:1:16],'FontSize',14,'FontName','Arial','YLim',[13,16]);
set(gca, 'XTick',[0.31:0.01:0.35],'FontSize',14,'FontName','Arial','XLim',[0.303,0.337]);
colormap(n4,nclCM(19));
caxis([1.3,1.9]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[1.3:0.3:1.9]);
set(get(h,'ylabel'),'string','Summer ER (PgC)','fontsize',14);
box on
% 为坐标区加标签
xlabel( 'Summer RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14);
ylabel( 'Summer TEM (°C)','FontName','Arial','FontSize',14);
zlabel( 'Summer ER (PgC)','FontName','Arial','FontSize',14);
grid on
view([-40.8648852383114 9.12669356141709]);
text('string','d','Units','normalized','position',[-0.16719463246739 1.03106688869709 0],'FontName','Arial','FontSize',14,'fontweight','bold')
plot3([0.303,0.303,0.337],[13,13,13],[1.3,1.9,1.9],'--','LineWidth',0.5,'color','k')
plot3([0.303,0.303],[13,16],[1.9,1.9],'--','LineWidth',0.5,'color','k')

d=axes('Position',[0.275070581981008 0.636489479512736 0.125788032787166 0.136212624584717]);
scatter(TEM_summer_list(:),ER_summer_list(:),'filled','MarkerFaceColor',[173,87,38]/255,'SizeData',60); hold on
scatter(TEM_summer_list(end),ER_summer_list(end),'filled','MarkerFaceColor',[0,107,172]/255,'SizeData',60); hold on
% 使用二项式拟合 (二次多项式)
[p,s] = polyfit(TEM_summer_list(:),ER_summer_list(:),2);  % p 是拟合系数，S 是拟合统计信息
x1=linspace(min(TEM_summer_list(:)),max(TEM_summer_list(:)));
% 计算拟合值
[y1, delta] = polyval(p, x1, s);  % y_fit 是拟合值，delta 是置信区间偏差
plot(x1,y1,'--','LineWidth',1.5,'color','k'); hold on
% 绘制置信区间阴影
fill([x1, fliplr(x1)], [y1 + delta, fliplr(y1 - delta)], ...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
R2=s.rsquared;
R2=roundn(R2,-2);
sig_text=['R^2 = ' num2str(R2)];
text('string',sig_text,'Units','normalized','position',[0.0989794543340237 0.128574645297961 0],'FontName','Arial','FontSize',10,'Color','k')
text(TEM_summer_list(end),ER_summer_list(end),'2023\rightarrow ','HorizontalAlignment','right','FontSize',10,'FontName','Arial','fontweight','bold','Color', [0,107,172]/255)
box on
set(gca, 'XTick',[13:1:16],'FontSize',10,'FontName','Arial','XLim',[13,16]);
set(gca, 'YTick',[1.2:0.3:1.8],'FontSize',10,'FontName','Arial','YLim',[1.2,1.9]);
xlabel( 'Summer TEM','FontName','Arial','FontSize',10);
ylabel( 'Summer ER','FontName','Arial','FontSize',10);


set(gcf,'unit','centimeters','position',[8.237361111111111,16.21366666666667,30.81513888888889,19.113500000000002]);
result=['E:\phd_file\Boreal_North_America\Result\V6\anomaly2023_summer_scatter_map_GCB2024.png']
% print(result,f,'-r600','-dpng');
