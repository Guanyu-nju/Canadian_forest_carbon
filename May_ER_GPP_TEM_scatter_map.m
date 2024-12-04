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
%% calculate mean value

for year=2015:2023

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_5.tif']);
    GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_5.tif']);
    ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;
    
    GCB_ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_5.tif']);
    GCB_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

    NEE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_' num2str(year) '_5.tif']);
    NEE_mean(:,:,year-2014)=NEE_mean_temp.*pixel_mask;

    GCB_NEE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_5.tif']);
    GCB_NEE_mean(:,:,year-2014)=GCB_NEE_mean_temp.*pixel_mask;
end
GPP_mean=nanmean(GPP_mean,3);
ER_mean=nanmean(ER_mean,3);
GCB_ER_mean=nanmean(GCB_ER_mean,3);
NEE_mean=nanmean(NEE_mean,3);
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);

GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_2023_5.tif").*pixel_mask;
ER2023=importdata("E:\phd_file\Boreal_North_America\ER\month\ER_2023_5.tif").*pixel_mask;
GCB_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_2023_5.tif").*pixel_mask;
NEE2023=importdata("E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_2023_5.tif").*pixel_mask;
GCB_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_2023_5.tif").*pixel_mask;

GPP_2023Anomaly=GPP2023-GPP_mean;
ER_2023Anomaly=ER2023-ER_mean;
GCB_ER_2023Anomaly=GCB_ER2023-GCB_ER_mean;
NEE_2023Anomaly=-NEE2023+NEE_mean;
GCB_NEE_2023Anomaly=-GCB_NEE2023+GCB_NEE_mean;

for year=2015:2023

    RZSM_mean_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_5.tif']);
    RZSM_may_mean(:,:,year-2014)=RZSM_mean_temp.*pixel_mask;

    TEM_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_5.tif']);
    TEM_may_mean(:,:,year-2014)=TEM_mean_temp.*pixel_mask;

end
RZSM_may_mean=nanmean(RZSM_may_mean,3);
TEM_may_mean=nanmean(TEM_may_mean,3);


RZSM_may_2023=importdata("E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_2023_5.tif").*pixel_mask;
TEM_may_2023=importdata("E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_2023_5.tif").*pixel_mask;

RZSM_2023_may_Anomaly=RZSM_may_2023-RZSM_may_mean;
TEM_2023_may_Anomaly=TEM_may_2023-TEM_may_mean;


for year=2015:2023

    RZSM_mean_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\spring\RZSM_spring_' num2str(year) '.tif']);
    RZSM_spring_mean(:,:,year-2014)=RZSM_mean_temp.*pixel_mask;

    TEM_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\spring\Temperature_spring_' num2str(year) '.tif']);
    TEM_spring_mean(:,:,year-2014)=TEM_mean_temp.*pixel_mask;

end
RZSM_spring_mean=nanmean(RZSM_spring_mean,3);
TEM_spring_mean=nanmean(TEM_spring_mean,3);


RZSM_spring_2023=importdata("E:\phd_file\Boreal_North_America\RZSM\spring\RZSM_spring_2023.tif").*pixel_mask;
TEM_spring_2023=importdata("E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\spring\Temperature_spring_2023.tif").*pixel_mask;

RZSM_2023_spring_Anomaly=RZSM_spring_2023-RZSM_spring_mean;
TEM_2023_spring_Anomaly=TEM_spring_2023-TEM_spring_mean;
%%
GPP=GPP_2023Anomaly(:);
ER=ER_2023Anomaly(:);
GCB_ER=GCB_ER_2023Anomaly(:);
NEE=NEE_2023Anomaly(:);
GCB_NEE=GCB_NEE_2023Anomaly(:);

RZSM_may=RZSM_2023_may_Anomaly(:);
TEM_may=TEM_2023_may_Anomaly(:);
RZSM_spring=RZSM_2023_spring_Anomaly(:);
TEM_spring=TEM_2023_spring_Anomaly(:);

RZSM_may(isnan(ER))=nan;
TEM_may(isnan(ER))=nan;
GPP(isnan(ER))=nan;
NEE(isnan(ER))=nan;
GCB_NEE(isnan(ER))=nan;
GCB_ER(isnan(ER))=nan;
RZSM_spring(isnan(ER))=nan;
TEM_spring(isnan(ER))=nan;

ER(isnan(ER))=[];
GCB_ER(isnan(GCB_ER))=[];
GCB_ER(isnan(GCB_ER))=[];

GCB_NEE(isnan(GCB_NEE))=[];
RZSM_may(isnan(RZSM_may))=[];
TEM_may(isnan(TEM_may))=[];
GPP(isnan(GPP))=[];
NEE(isnan(NEE))=[];
RZSM_spring(isnan(RZSM_spring))=[];
TEM_spring(isnan(TEM_spring))=[];

f=figure
t = tiledlayout(6,5);

t.TileSpacing = 'compact';
t.Padding = 'compact';

% part 1 *********************************************************************************
n1=nexttile([2,3])
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GPP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n1,nclCM(399));
caxis([-40,40]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('May GPP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.0597839402568452 1.08038697114719 0],'FontName','Arial','FontSize',14,'fontweight','bold')

n5=nexttile([2,2])
scatter(TEM_spring,GPP,[],RZSM_spring,'filled','SizeData',80); hold on
colormap(n5,nclCM(226));
caxis([-0.06,0.06]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.06:0.03:0.06]);
set(get(h,'ylabel'),'string','Spring RZSM anomaly (m^{3} m^{-3})','fontsize',14);
box on
ylabel('May GPP anomalies (gC m^{-2})','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Spring TEM anomalies (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[-2:2:6],'FontSize',14,'FontName','Arial','XLim',[-2,6]);
% set(gca, 'YTick',[-40:20:40],'FontSize',14,'FontName','Arial','YLim',[-40,40]);
text('string','b','Units','normalized','position',[-0.353103150572238 1.04917033415992 0],'FontName','Arial','FontSize',14,'fontweight','bold')
[r,p]=corrcoef(TEM_spring,GPP);
p=p(1,2);
r=r(1,2);
sig_text=['R \rm= ' num2str(roundn(r,-2))];
if p<0.01
    sig_text=[sig_text,'**'];
elseif p<0.05 & p>0.01
    sig_text=[sig_text,'*'];
end
text('string',sig_text,'Units','normalized','position',[0.692861929523356 0.0586797598440925 0],'FontName','Arial','FontSize',14)

% part 2 *********************************************************************************
n2=nexttile([2,3])
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n2,flipud(nclCM(399)));
caxis([-40,40]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('May ER anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.0597839402568452 1.08038697114719 0],'FontName','Arial','FontSize',14,'fontweight','bold')

n6=nexttile([2,2])
scatter(TEM_spring,ER,[],RZSM_spring,'filled','SizeData',80); hold on
colormap(n6,nclCM(226));
caxis([-0.06,0.06]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.06:0.03:0.06]);
set(get(h,'ylabel'),'string','Spring RZSM anomaly (m^{3} m^{-3})','fontsize',14);
box on
ylabel('May ER anomalies (gC m^{-2})','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Spring TEM anomalies (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[-2:2:6],'FontSize',14,'FontName','Arial','XLim',[-2,6]);
% set(gca, 'YTick',[-40:20:40],'FontSize',14,'FontName','Arial','YLim',[-40,40]);
[r,p]=corrcoef(TEM_spring,ER);
p=p(1,2);
r=r(1,2);
sig_text=['R \rm= ' num2str(roundn(r,-2))];
if p<0.01
    sig_text=[sig_text,'**'];
elseif p<0.05 & p>0.01
    sig_text=[sig_text,'*'];
end
text('string',sig_text,'Units','normalized','position',[0.692861929523356 0.0586797598440925 0],'FontName','Arial','FontSize',14)
text('string','d','Units','normalized','position',[-0.353103150572238 1.04917033415992 0],'FontName','Arial','FontSize',14,'fontweight','bold')
text('string','GCASv2','Units','normalized','position',[0.688563935254015 0.942580679545433 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 3 *********************************************************************************
n2=nexttile([2,3])
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n2,flipud(nclCM(399)));
caxis([-40,40]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('May ER anomalies (GCB2024)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','e','Units','normalized','position',[-0.0597839402568452 1.08038697114719 0],'FontName','Arial','FontSize',14,'fontweight','bold')

n6=nexttile([2,2])
scatter(TEM_spring,GCB_ER,[],RZSM_spring,'filled','SizeData',80); hold on
colormap(n6,nclCM(226));
caxis([-0.06,0.06]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.06:0.03:0.06]);
set(get(h,'ylabel'),'string','Spring RZSM anomaly (m^{3} m^{-3})','fontsize',14);
box on
ylabel('May ER anomalies (gC m^{-2})','FontName','Arial','FontSize',14)  % 设置x轴标题
xlabel('Spring TEM anomalies (°C)','FontName','Arial','FontSize',14)  % 设置x轴标题
set(gca, 'XTick',[-2:2:6],'FontSize',14,'FontName','Arial','XLim',[-2,6]);
set(gca, 'YTick',[-40:20:40],'FontSize',14,'FontName','Arial','YLim',[-40,40]);
[r,p]=corrcoef(TEM_spring,GCB_ER);
p=p(1,2);
r=r(1,2);
sig_text=['R \rm= ' num2str(roundn(r,-2))];
if p<0.01
    sig_text=[sig_text,'**'];
elseif p<0.05 & p>0.01
    sig_text=[sig_text,'*'];
end
text('string',sig_text,'Units','normalized','position',[0.692861929523356 0.0586797598440925 0],'FontName','Arial','FontSize',14)
text('string','f','Units','normalized','position',[-0.353103150572238 1.04917033415992 0],'FontName','Arial','FontSize',14,'fontweight','bold')
text('string','GCB2024','Units','normalized','position',[0.688563935254015 0.942580679545433 0],'FontName','Arial','FontSize',14,'fontweight','bold')

set(gcf,'unit','centimeters','position',[26.431875000000005,6.1595,30.532916666666672,27.574875000000006]);
result=['E:\phd_file\Boreal_North_America\Result\V6\May_ER_GPP_TEM_scatter_map.png']
% print(result,f,'-r600','-dpng');
