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

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_6.tif']);
    June_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_6.tif']);
    June_GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_6.tif']);
    June_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_7.tif']);
    July_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_7.tif']);
    July_GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_7.tif']);
    July_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_8.tif']);
    Aug_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_8.tif']);
    Aug_GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_8.tif']);
    Aug_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;
end
June_NEE_mean=nanmean(June_NEE_mean,3);
June_GPP_mean=nanmean(June_GPP_mean,3);
June_ER_mean=nanmean(June_ER_mean,3);

July_NEE_mean=nanmean(July_NEE_mean,3);
July_GPP_mean=nanmean(July_GPP_mean,3);
July_ER_mean=nanmean(July_ER_mean,3);

Aug_NEE_mean=nanmean(Aug_NEE_mean,3);
Aug_GPP_mean=nanmean(Aug_GPP_mean,3);
Aug_ER_mean=nanmean(Aug_ER_mean,3);

June_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_2023_6.tif").*pixel_mask;
June_GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_2023_6.tif").*pixel_mask;
June_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_2023_6.tif").*pixel_mask;

July_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_2023_7.tif").*pixel_mask;
July_GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_2023_7.tif").*pixel_mask;
July_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_2023_7.tif").*pixel_mask;

Aug_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_2023_8.tif").*pixel_mask;
Aug_GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_2023_8.tif").*pixel_mask;
Aug_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_2023_8.tif").*pixel_mask;

June_NEP_2023Anomaly=-June_NEE2023+June_NEE_mean;
June_GPP_2023Anomaly=June_GPP2023-June_GPP_mean;
June_ER_2023Anomaly=June_ER2023-June_ER_mean;

July_NEP_2023Anomaly=-July_NEE2023+July_NEE_mean;
July_GPP_2023Anomaly=July_GPP2023-July_GPP_mean;
July_ER_2023Anomaly=July_ER2023-July_ER_mean;

Aug_NEP_2023Anomaly=-Aug_NEE2023+Aug_NEE_mean;
Aug_GPP_2023Anomaly=Aug_GPP2023-Aug_GPP_mean;
Aug_ER_2023Anomaly=Aug_ER2023-Aug_ER_mean;

%%
f=figure
t = tiledlayout(2,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';


% part 1 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,June_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('June NEP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 2 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,July_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('July NEP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 3 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,Aug_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('August NEP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 7 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,June_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,flipud(nclCM(399)));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('June ER anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 8 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,July_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,flipud(nclCM(399)));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('July ER anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','e','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 9 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,Aug_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,flipud(nclCM(399)));

caxis([-80,80]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-80:40:80]);
set(get(h,'ylabel'),'string','gC m^{-2} month^{-1}','fontsize',14);
title('August ER anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','f','Units','normalized','position',[-0.0697588532029519 1.12202323328086 0],'FontName','Arial','FontSize',14,'fontweight','bold')

set(gcf,'unit','centimeters','position',[10.054166666666669,13.096875000000002,45.82583333333334,19.54212499999999]);
result=['E:\phd_file\Boreal_North_America\Result\V6\Flux_anomaly2023_June_July_Aug_map_GCB2024.png']
% print(result,f,'-r600','-dpng');
