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

for year=2015:2023

    Con_April=importdata(['E:\phd_file\Boreal_North_America\OCO2-XCO2\OCO2-XCO2_' num2str(year) '_4.tif']);
    April_mask=Con_April;
    April_mask(~isnan(April_mask))=1;
    Con_August=importdata(['E:\phd_file\Boreal_North_America\OCO2-XCO2\OCO2-XCO2_' num2str(year) '_8.tif']);

    April_mask(isnan(Con_August))=nan;
    Con_diff=nansum(nansum(Con_April.*pixel_mask.*area_grid.*April_mask))/(nansum(nansum(area_grid.*April_mask)))-nansum(nansum(Con_August.*pixel_mask.*area_grid.*April_mask))/nansum(nansum(area_grid.*April_mask));


    Data_list(year-2014)=Con_diff;


end

f=figure
t = tiledlayout(1,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
b1=bar(Data_list,'FaceColor','flat');hold on
for i=1:8
    b1.CData(i,:) = [55,140,230]/255;
end
b1.CData(9,:) = [237,145,29]/255;
ylim([6,10.5])
set(gca,'YTick', [6:1:10]);
ylabel('ppm','FontName','Arial','FontSize',14);
xlim([0.25,9.75])
set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2015:1:2023],'FontName','Arial','fontsize',14)
xlabel('Year','FontName','Arial','FontSize',14)
title('Difference in XCO_2 concentration between April and August','FontName','Arial','FontSize',14,'fontweight','bold')

set(gcf,'unit','centimeters','position',[25.26770833333334,15.853833333333334,16.91745833333333,9.731375000000007]);
result=['E:\phd_file\Boreal_North_America\Result\V6\CO2_diff_bar.png']
% print(result,f,'-r600','-dpng');