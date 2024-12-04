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

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\year\NEE_' num2str(year) '.tif']);
    NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\year\ER_' num2str(year) '.tif']);
    ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
NEE_mean=nanmean(NEE_mean,3);
GPP_mean=nanmean(GPP_mean,3);
ER_mean=nanmean(ER_mean,3);



NEE2023=importdata("E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\year\NEE_2023.tif").*pixel_mask;
GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_2023.tif").*pixel_mask;
ER2023=importdata("E:\phd_file\Boreal_North_America\ER\year\ER_2023.tif").*pixel_mask;


NEP_2023Anomaly=-NEE2023+NEE_mean;
GPP_2023Anomaly=GPP2023-GPP_mean;
ER_2023Anomaly=ER2023-ER_mean;
%% calculate mean value (GCB2024)

for year=2015:2023

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\year\NEE_' num2str(year) '.tif']);
    GCB_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GCB_GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\year\ER_' num2str(year) '.tif']);
    GCB_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);
GCB_GPP_mean=nanmean(GCB_GPP_mean,3);
GCB_ER_mean=nanmean(GCB_ER_mean,3);



GCB_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\year\NEE_2023.tif").*pixel_mask;
GCB_GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_2023.tif").*pixel_mask;
GCB_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\year\ER_2023.tif").*pixel_mask;


GCB_NEP_2023Anomaly=-GCB_NEE2023+GCB_NEE_mean;
GPP_2023Anomaly=GCB_GPP2023-GCB_GPP_mean;
GCB_ER_2023Anomaly=GCB_ER2023-GCB_ER_mean;
%%
% 面积占比
count_sum=sum(sum(~isnan(GCB_NEP_2023Anomaly)));
count1=sum(sum(GCB_NEP_2023Anomaly>0));
result=count1/(count_sum)
%% GCAS
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;

% 时间序列均值
for month=1:12


    for year=2015:2023
        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));

        BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);
        NEE_std_list(year-2014,month,1)=nansum(nansum(BEPS_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
        NEE_std_list(year-2014,month,2)=nansum(nansum(BEPS_GFED_NEE.*pixel_mask.*area_grid/(10^15)));
        NEE_std_list(year-2014,month,3)=nansum(nansum(CASA_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
        NEE_std_list(year-2014,month,4)=nansum(nansum(CASA_GFED_NEE.*pixel_mask.*area_grid/(10^15)));

        GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
        GPP_list(year-2014,month)=nansum(nansum(GPP_mean_temp.*pixel_mask.*area_grid/(10^15)));

        GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
        Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

        GPP_std_list(year-2014,month,1)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
        GPP_std_list(year-2014,month,2)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));


        ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_' num2str(month) '.tif']);
        ER_list(year-2014,month)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));


    end

end

% 2023年数据
for month=1:12


    year=2023;
    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    NEE_2023(month)=nansum(nansum(NNE_2023_temp.*pixel_mask.*area_grid/(10^15)));

    BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);
    NEE_2023_list(1,month)=nansum(nansum(BEPS_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_2023_list(2,month)=nansum(nansum(BEPS_GFED_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_2023_list(3,month)=nansum(nansum(CASA_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_2023_list(4,month)=nansum(nansum(CASA_GFED_NEE.*pixel_mask.*area_grid/(10^15)));


    GPP_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
    GPP_2023(month)=nansum(nansum(GPP_2023_temp.*pixel_mask.*area_grid/(10^15)));

    GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
        Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

        GPP_2023_list(1,month)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
        GPP_2023_list(2,month)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));


    ER_2023_temp=importdata(['E:\phd_file\Boreal_North_America\ER\month\ER_' num2str(year) '_' num2str(month) '.tif']);
    ER_2023(month)=nansum(nansum(ER_2023_temp.*pixel_mask.*area_grid/(10^15)));


end

% NEP距平均值与标准差
NEP_anomaly_mean=-NEE_2023+nanmean(NEE_list);
for i=1:size(NEE_std_list,3)

    temp_list=NEE_std_list(:,:,i);
    temp2023=NEE_2023_list(i,:);
    NEP_anomaly_std(i,:)=-temp2023+nanmean(temp_list);

end
NEP_anomaly_std=nanstd(NEP_anomaly_std);

% GPP距平均值
GPP_anomaly_mean=GPP_2023-nanmean(GPP_list);
for i=1:size(GPP_std_list,3)

    temp_list=GPP_std_list(:,:,i);
    temp2023=GPP_2023_list(i,:);
    GPP_anomaly_std(i,:)=temp2023-nanmean(temp_list);

end
% GPP_anomaly_std=nanstd(GPP_anomaly_std);
for i=1:length(GPP_anomaly_std)
    GPP_temp=GPP_std_list(:,i);
    GPP_std(i)=std(bootstrp(1000,@mean,GPP_temp));
end
GPP_anomaly_std=GPP_std;

% ER距平均值与标准差
ER_anomaly_mean=ER_2023-nanmean(ER_list);

ER_anomaly_std=sqrt(power(NEP_anomaly_std,2)+power(GPP_anomaly_std,2));


%% GCB2024

% 时间序列均值
for month=1:12


    for year=2015:2023
        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GCB_NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));


        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);


        GCB_NEE_std_list(year-2014,month,1)=nansum(nansum(A_NEE.*pixel_mask.*area_grid/(10^15)));
        GCB_NEE_std_list(year-2014,month,2)=nansum(nansum(B_NEE.*pixel_mask.*area_grid/(10^15)));
        GCB_NEE_std_list(year-2014,month,3)=nansum(nansum(C_NEE.*pixel_mask.*area_grid/(10^15)));
        GCB_NEE_std_list(year-2014,month,4)=nansum(nansum(D_NEE.*pixel_mask.*area_grid/(10^15)));
   
        ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_' num2str(month) '.tif']);
        GCB_ER_list(year-2014,month)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    end

end

% 2023年数据
for month=1:12


    year=2023;
    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    GCB_NEE_2023(month)=nansum(nansum(NNE_2023_temp.*pixel_mask.*area_grid/(10^15)));


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);


    GCB_NEE_2023_list(1,month)=nansum(nansum(A_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_2023_list(2,month)=nansum(nansum(B_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_2023_list(3,month)=nansum(nansum(C_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_2023_list(4,month)=nansum(nansum(D_NEE.*pixel_mask.*area_grid/(10^15)));

    ER_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_' num2str(month) '.tif']);
    GCB_ER_2023(month)=nansum(nansum(ER_2023_temp.*pixel_mask.*area_grid/(10^15)));


end

% NEP距平均值与标准差
GCB_NEP_anomaly_mean=-GCB_NEE_2023+nanmean(GCB_NEE_list);
for i=1:size(GCB_NEE_std_list,3)

    temp_list=GCB_NEE_std_list(:,:,i);
    temp2023=GCB_NEE_2023_list(i,:);
    GCB_NEP_anomaly_std(i,:)=-temp2023+nanmean(temp_list);

end
GCB_NEP_anomaly_std=nanstd(GCB_NEP_anomaly_std);

% ER距平均值与标准差
GCB_ER_anomaly_mean=GCB_ER_2023-nanmean(GCB_ER_list);
GCB_ER_anomaly_std=sqrt(power(GCB_NEP_anomaly_std,2)+power(GPP_anomaly_std,2));

%%
Fire2023=importdata("E:\phd_file\Boreal_North_America\fire emission\mean\year\Fire_2023.tif").*pixel_mask;
for year=2015:2023

    Fire_mean_temp=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\year\Fire_' num2str(year) '.tif']);
    Fire_mean(:,:,year-2014)=Fire_mean_temp.*pixel_mask;


end
Fire_mean=nanmean(Fire_mean,3);
Fire_2023Anomaly=Fire2023-Fire_mean;

%%
f=figure
t = tiledlayout(3,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';


% part 1 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 2 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (GCB2024)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')


% part 3 *********************************************************************************
n4=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GPP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n4,nclCM(399));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('GPP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 4 *********************************************************************************
n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n6,flipud(nclCM(399)));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('ER anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 5 *********************************************************************************
n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n6,flipud(nclCM(399)));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('ER anomalies (GCB2024)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','e','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 6 *********************************************************************************
n7=nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,Fire_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(n7,nclCM(432));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('Fire CO_2 emission anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','f','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',14,'fontweight','bold')


% part 7 ************************************************************************<0*********
nexttile
x=1:length(NEP_anomaly_mean);
plot([0,20],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[157,157,157]/255);hold on
p1=plot(NEP_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);
fill([x, fliplr(x)], [NEP_anomaly_mean, fliplr(NEP_anomaly_mean+NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [NEP_anomaly_mean, fliplr(NEP_anomaly_mean-NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p2=plot(GPP_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[17,119,51]/255);
fill([x, fliplr(x)], [GPP_anomaly_mean, fliplr(GPP_anomaly_mean+GPP_anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [GPP_anomaly_mean, fliplr(GPP_anomaly_mean-GPP_anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p3=plot(ER_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[193,0,1]/255);hold on
fill([x, fliplr(x)], [ER_anomaly_mean, fliplr(ER_anomaly_mean+ER_anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [ER_anomaly_mean, fliplr(ER_anomaly_mean-ER_anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);


ylim([-0.17,0.21])
set(gca,'YTick', [-0.1:0.1:0.2]);
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','g','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','ER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
title('Monthly flux anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')

% part 8 *********************************************************************************
nexttile
x=1:length(GCB_NEP_anomaly_mean);
plot([0,20],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[157,157,157]/255);hold on
p1=plot(GCB_NEP_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);
fill([x, fliplr(x)], [GCB_NEP_anomaly_mean, fliplr(GCB_NEP_anomaly_mean+GCB_NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [GCB_NEP_anomaly_mean, fliplr(GCB_NEP_anomaly_mean-GCB_NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p2=plot(GPP_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[17,119,51]/255);
fill([x, fliplr(x)], [GPP_anomaly_mean, fliplr(GPP_anomaly_mean+GPP_anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [GPP_anomaly_mean, fliplr(GPP_anomaly_mean-GPP_anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p3=plot(GCB_ER_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[193,0,1]/255);hold on
fill([x, fliplr(x)], [GCB_ER_anomaly_mean, fliplr(GCB_ER_anomaly_mean+GCB_ER_anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [GCB_ER_anomaly_mean, fliplr(GCB_ER_anomaly_mean-GCB_ER_anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);

ylim([-0.17,0.21])
set(gca,'YTick', [-0.1:0.1:0.2]);
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','h','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','ER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
title('Monthly flux anomalies (GCB2024)','FontName','Arial','FontSize',14,'fontweight','bold')


set(gcf,'unit','centimeters','position',[10.054166666666667,8.128,45.825833333333335,24.510999999999996]);
result=['E:\phd_file\Boreal_North_America\Result\V6\Flux_Aomalies_year_map.png']
% print(result,f,'-r600','-dpng');
