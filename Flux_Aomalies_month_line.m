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
for month=1:12


    for year=2015:2023

        snow_melt_temp=importdata(['E:\phd_file\Boreal_North_America\snowmelt\month\ERA5_snowmelt\global\snow_melt_' num2str(year) '_' num2str(month) '.tif']);
        snow_melt_list(year-2014,month)=nansum(nansum(snow_melt_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        Temperature_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_' num2str(month) '.tif']);
        Temperature_list(year-2014,month)=nansum(nansum(Temperature_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        PPT_temp=importdata(['E:\phd_file\Boreal_North_America\PPT\month\ERA5_PPT\globe\ERA5_PPT_' num2str(year) '_' num2str(month) '.tif']);
        PPT_list(year-2014,month)=nansum(nansum(PPT_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_' num2str(month) '.tif']);
        RZSM_list(year-2014,month)=nansum(nansum(RZSM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        VPD_temp=importdata(['E:\phd_file\Boreal_North_America\VPD\VPD\month\globe\VPD_' num2str(year) '_' num2str(month) '.tif']);
        VPD_list(year-2014,month)=nansum(nansum(VPD_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        SR_temp=importdata(['E:\phd_file\Boreal_North_America\Solar_radiation\month\globe\PAR_' num2str(year) '_' num2str(month) '.tif']);
        SR_list(year-2014,month)=nansum(nansum(SR_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));


    end

end

% 时间序列均值
for month=1:12

    year=2023;

    snow_melt_temp=importdata(['E:\phd_file\Boreal_North_America\snowmelt\month\ERA5_snowmelt\global\snow_melt_' num2str(year) '_' num2str(month) '.tif']);
    snow_melt_2023(month)=nansum(nansum(snow_melt_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    Temperature_temp=importdata(['E:\phd_file\Boreal_North_America\ERA5_land_air_temperature\month\month_globe\ERA5_land_air_temperature_' num2str(year) '_' num2str(month) '.tif']);
    Temperature_2023(month)=nansum(nansum(Temperature_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    PPT_temp=importdata(['E:\phd_file\Boreal_North_America\PPT\month\ERA5_PPT\globe\ERA5_PPT_' num2str(year) '_' num2str(month) '.tif']);
    PPT_2023(month)=nansum(nansum(PPT_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    RZSM_temp=importdata(['E:\phd_file\Boreal_North_America\RZSM\month\globe\ERA5_RZSM_' num2str(year) '_' num2str(month) '.tif']);
    RZSM_2023(month)=nansum(nansum(RZSM_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

    VPD_temp=importdata(['E:\phd_file\Boreal_North_America\VPD\VPD\month\globe\VPD_' num2str(year) '_' num2str(month) '.tif']);
    VPD_2023(month)=nansum(nansum(VPD_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

        SR_temp=importdata(['E:\phd_file\Boreal_North_America\Solar_radiation\month\globe\PAR_' num2str(year) '_' num2str(month) '.tif']);
    SR_2023(month)=nansum(nansum(SR_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));

end


%%
f=figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% part 1 *********************************************************************************
nexttile
p1=plot(nanmean(Temperature_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(Temperature_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
ylim([-20,20])
set(gca,'YTick', [-20:10:25]);
ylabel('TEM (°C)','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','a','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2023'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

% part 2 *********************************************************************************
nexttile
p1=plot(nanmean(PPT_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(PPT_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
ylim([28,100])
set(gca,'YTick', [20:20:120]);
ylabel('PRE (mm)','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','b','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2022'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

% part 1 *********************************************************************************
nexttile
p1=plot(nanmean(RZSM_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(RZSM_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
ylim([0.28,0.38])
set(gca,'YTick', [0.28:0.02:0.38]);
ylabel('RZSM (m^{3} m^{-3})','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','c','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2023'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

% part 1 *********************************************************************************
nexttile
p1=plot(nanmean(snow_melt_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(snow_melt_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
ylim([0,0.16])
set(gca,'YTick', [0:0.04:0.16]);
ylabel('Snowmelt (m of water equivalent)','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','d','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2023'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

% part 5 *********************************************************************************
nexttile
p1=plot(nanmean(VPD_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(VPD_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
% ylim([-1.5,1])
% set(gca,'YTick', [-1.5:0.5:1]);
ylabel('VPD (kPa)','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','e','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2023'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

% part 6 *********************************************************************************
nexttile
p1=plot(nanmean(SR_list),'-','LineWidth',2,'MarkerSize',15,'color','k'); hold on
p2=plot(SR_2023,'.-','LineWidth',2,'MarkerSize',20,'color','r');
ylim([0,700])
set(gca,'YTick', [0:200:600]);
ylabel('SR (MJ m^-^2)','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','f','Units','normalized','position',[-0.121211776822882 1.06809391674095 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 p2],{'2015-2023','2023'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)

set(gcf,'unit','centimeters','position',[26.431875000000005,2.899833333333333,34.44345833333333,30.83454166666668]);
result=['E:\phd_file\Boreal_North_America\Result\V6\Flux_Aomalies_month_line.png']
% print(result,f,'-r600','-dpng');



