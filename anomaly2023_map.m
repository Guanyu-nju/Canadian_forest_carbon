%% Mask create
clc
clear
% Canada mask
pixel_mask=importdata("E:\phd_file\Boreal_North_America\Canada_mask_0.1degree.tif");
land_mask=importdata("E:\phd_file\Boreal_North_America\Global_PRE_TEM\TEM\ERA5_land_air_temperature\ERA5_land_air_temperature_2023.tif");
land_mask(~isnan(land_mask))=1;
pixel_mask=pixel_mask.*land_mask;

% 创建环境
data=geotiffread("E:\phd_file\Boreal_North_America\Canada_mask_0.1degree.tif");
info=geotiffinfo("E:\phd_file\Boreal_North_America\Canada_mask_0.1degree.tif");
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

clearvars -except lon1 lat1 land_mask pixel_mask bou_canX bou_canY m n
%% 2023 zscore

Temperature_year2023=importdata("E:\phd_file\Boreal_North_America\Global_PRE_TEM\TEM\ERA5_land_air_temperature\ERA5_land_air_temperature_2023.tif");
PPT_year2023=importdata("E:\phd_file\Boreal_North_America\Global_PRE_TEM\PRE\ERA5_land_PPT\ERA5_land_PPT_2023.tif");


for year=2000:2023

    TEM_year_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\TEM\ERA5_land_air_temperature\ERA5_land_air_temperature_' num2str(year) '.tif']);
    TEM_year_mean(:,:,year-1999)=TEM_year_temp;

    PPT_year_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\PRE\ERA5_land_PPT\ERA5_land_PPT_' num2str(year) '.tif']);
    PPT_year_mean(:,:,year-1999)=PPT_year_temp;

end
[~,TEM_sort]=sort(TEM_year_mean,3);
TEM_year_sort=TEM_sort(:,:,end);
[~,PPT_sort]=sort(PPT_year_mean,3);
PPT_year_sort=PPT_sort(:,:,end);

TEM_year_mean=nanmean(TEM_year_mean,3);
TEM_absAnomaly_year=Temperature_year2023-TEM_year_mean;
PPT_year_mean=nanmean(PPT_year_mean,3);
PPT_absAnomaly_year=(PPT_year2023-PPT_year_mean)./PPT_year_mean*100;

%% calculate mean value

area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter_0.1degree.tif")*1000000.*land_mask;
Canada_area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter_0.1degree.tif")*1000000.*pixel_mask;


% 时间序列均值

for year=2000:2023

    Temperature_year_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\TEM\ERA5_land_air_temperature\ERA5_land_air_temperature_' num2str(year) '.tif']);
    Temperature_year_list(year-1999)=nansum(nansum(Temperature_year_temp.*area_grid))/(nansum(nansum(area_grid)));

    Temperature_canada_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\TEM\ERA5_land_air_temperature\ERA5_land_air_temperature_' num2str(year) '.tif']);
    Temperature_canada_list(year-1999)=nansum(nansum(Temperature_canada_temp.*Canada_area_grid))/(nansum(nansum(Canada_area_grid)));

    PPT_year_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\PRE\ERA5_land_PPT\ERA5_land_PPT_' num2str(year) '.tif']);
    PPT_year_list(year-1999)=nansum(nansum(PPT_year_temp.*area_grid))/(nansum(nansum(area_grid)));

    PPT_canada_temp=importdata(['E:\phd_file\Boreal_North_America\Global_PRE_TEM\PRE\ERA5_land_PPT\ERA5_land_PPT_' num2str(year) '.tif']);
    PPT_canada_list(year-1999)=nansum(nansum(PPT_canada_temp.*Canada_area_grid))/(nansum(nansum(Canada_area_grid)));

end

Temperature_year_list=Temperature_year_list-nanmean(Temperature_year_list);
Temperature_canada_list=Temperature_canada_list-nanmean(Temperature_canada_list);
PPT_year_list=(PPT_year_list-nanmean(PPT_year_list))./nanmean(PPT_year_list)*100;
PPT_canada_list=(PPT_canada_list-nanmean(PPT_canada_list))./nanmean(PPT_canada_list)*100;

%%

f=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% part 1 *********************************************************************************
n1=nexttile
m_proj('Equidistant','long',[-180 180],'lat',[-60 90]);
m_pcolor(lon1,lat1,TEM_absAnomaly_year,'linestyle','none');
m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-180:60:180],'ytick',[-30:30:90],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-180:60:180],'yticklabels',[-30:30:90],'linestyle','none'); %对齐格网
colormap(n1,nclCM(79));
caxis([-3,3]);
h=colorbar('location','southoutside','FontName','Arial','FontSize',14,'ytick',[-3:1:3]);
set(get(h,'ylabel'),'string','Absolute anomalies (°C)','fontsize',14);
title('2023 TEM anomalies relative to 2000-2023','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.12721755945316 1.09694973692137 0],'FontName','Arial','FontSize',14,'fontweight','bold')

nexttile
b1=bar(Temperature_year_list,'FaceColor',[217,222,231]/255);hold on
p1=plot(Temperature_canada_list,'-','LineWidth',2,'MarkerSize',15,'color',[226,85,45]/255);hold on
ylim([-1.5,2.7])
set(gca,'YTick', [-1:1:2]);
ylabel('TEM anomalies (°C)','FontName','Arial','FontSize',14);
% xlim([0.5,24.5])
set(gca,'XTick',[1:5:21,24],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2000:5:2020,2023],'FontName','Arial','fontsize',14)
text('string','b','Units','normalized','position',[-0.100960171626526 1.06422561551247 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1,b1],{'Canada','Global'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
% text('string','RZSM','Units','normalized','position',[0.784819902167739 0.908688370089802 0],'FontName','Arial','FontSize',14,'fontweight','bold')
xlabel('Year','FontName','Arial','FontSize',14)

% part 2 *********************************************************************************
n1=nexttile
m_proj('Equidistant','long',[-180 180],'lat',[-60 90]);
m_pcolor(lon1,lat1,PPT_absAnomaly_year,'linestyle','none');
m_line(bou_canX,bou_canY,'linewidth',1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-180:60:180],'ytick',[-30:30:90],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-180:60:180],'yticklabels',[-30:30:90],'linestyle','none'); %对齐格网
colormap(n1,flipud(nclCM(79)));
caxis([-40,40]);
h=colorbar('location','southoutside','FontName','Arial','FontSize',14,'ytick',[-40:10:40]);
set(get(h,'ylabel'),'string','Relative anomalies (%)','fontsize',14);
title('2023 PRE anomalies relative to 2000-2023','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.12721755945316 1.09694973692137 0],'FontName','Arial','FontSize',14,'fontweight','bold')

nexttile
b1=bar(PPT_year_list,'FaceColor',[217,222,231]/255);hold on
p1=plot(PPT_canada_list,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
ylim([-10,15])
set(gca,'YTick', [-10:5:15]);
ylabel('PRE anomalies (%)','FontName','Arial','FontSize',14);
% xlim([0.5,24.5])
set(gca,'XTick',[1:5:21,24],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2000:5:2020,2023],'FontName','Arial','fontsize',14)
text('string','d','Units','normalized','position',[-0.100960171626526 1.06422561551247 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1,b1],{'Canada','Global'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
% text('string','RZSM','Units','normalized','position',[0.784819902167739 0.908688370089802 0],'FontName','Arial','FontSize',14,'fontweight','bold')
xlabel('Year','FontName','Arial','FontSize',14)

set(gcf,'unit','centimeters','position',[19.658541666666668,4.614333333333334,34.655125,20.494625]);
result=['E:\phd_file\Boreal_North_America\Result\V6\anomaly2023_map']
% print(result,f,'-r600','-dpng');
