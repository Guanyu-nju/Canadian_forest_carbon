clc
clear
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

%%
etl_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_etl_monthly_co2_SG.xlsx");
% 提取每年4月和8月的数据
data_april = etl_data(etl_data.months == 4, :);
data_august = etl_data(etl_data.months == 8, :);
co2_diff=data_april{:,end}-data_august{:,end};
co2_diff=[co2_diff(1:5);nan;co2_diff(6:end)];
etl_data{:,end}=detrend(etl_data{:,end},"omitnan");
% 计算每个年份有多少个月的数据
months_per_year = varfun(@numel, etl_data, 'InputVariables', 'months', 'GroupingVariables', 'years');
% 仅保留有12个月数据的年份
full_years = months_per_year.years(months_per_year.numel_months == 12);
filtered_table = etl_data(ismember(etl_data.years, full_years), :);
% 提取2023年数据
data_2023 = filtered_table(filtered_table.years == 2023 & filtered_table.months >= 3 & filtered_table.months <= 11, :);
% 提取其他年份数据（剔除2023年）
% other_years_data = filtered_table(filtered_table.years ~= 2023 & filtered_table.months >= 3 & filtered_table.months <= 11, :);
other_years_data = filtered_table( filtered_table.months >= 3 & filtered_table.months <= 11, :);

% 计算其他年份的每个月的均值和标准差
monthly_avg = varfun(@mean, other_years_data, 'InputVariables', 'mean_filtered_concentration', ...
    'GroupingVariables', 'months');
monthly_std = varfun(@std, other_years_data, 'InputVariables', 'mean_filtered_concentration', ...
    'GroupingVariables', 'months');
co2_2023=data_2023{:,end};
co2_mean=monthly_avg{:,end}';
co2_std=monthly_std{:,end}';
co2_anomaly=[nanmean(co2_2023(1:3))-nanmean(co2_mean(1:3)),nanmean(co2_2023(4:6))-nanmean(co2_mean(4:6)),nanmean(co2_2023(7:9))-nanmean(co2_mean(7:9))];

% 计算pixel对应数值
lat=54.9538;
lon=-112.4666;
[~,min_lat_index] = min(abs(lat_-lat));
[~,min_lon_index] = min(abs(lon_-lon));
% 时间序列均值
for month=1:12

    for year=2015:2023
        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        NEE_list(year-2014,month)=NNE_mean_temp(min_lat_index,min_lon_index);

        BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);
        NEE_std_list(year-2014,month,1)=BEPS_GFAS_NEE(min_lat_index,min_lon_index);
        NEE_std_list(year-2014,month,2)=BEPS_GFED_NEE(min_lat_index,min_lon_index);
        NEE_std_list(year-2014,month,3)=CASA_GFAS_NEE(min_lat_index,min_lon_index);
        NEE_std_list(year-2014,month,4)=CASA_GFED_NEE(min_lat_index,min_lon_index);

        Fire_mean_temp=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Fire_list(year-2014,month)=Fire_mean_temp(min_lat_index,min_lon_index);

    end

end
% 2023年数据
for month=1:12


    year=2023;
    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    NEE_2023(month)=NNE_2023_temp(min_lat_index,min_lon_index);

    BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);
    NEE_2023_list(1,month)=BEPS_GFAS_NEE(min_lat_index,min_lon_index);
    NEE_2023_list(2,month)=BEPS_GFED_NEE(min_lat_index,min_lon_index);
    NEE_2023_list(3,month)=CASA_GFAS_NEE(min_lat_index,min_lon_index);
    NEE_2023_list(4,month)=CASA_GFED_NEE(min_lat_index,min_lon_index);

    Fire_2023_temp=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
    Fire_2023(month)=Fire_2023_temp(min_lat_index,min_lon_index);


end
% NEP距平均值与标准差
NEP_anomaly_mean=-NEE_2023+nanmean(NEE_list);
NEP_anomaly_mean=NEP_anomaly_mean(3:11);
for i=1:size(NEE_std_list,3)

    temp_list=NEE_std_list(:,:,i);
    temp2023=NEE_2023_list(i,:);
    NEP_anomaly_std(i,:)=-temp2023+nanmean(temp_list);

end
NEP_anomaly_std=nanstd(NEP_anomaly_std);
NEP_anomaly_std=NEP_anomaly_std(3:11);
% ER距平均值与标准差
Fire_anomaly_mean=Fire_2023-nanmean(Fire_list);



%%
f=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile([2,4])
b1=bar(co2_diff,'FaceColor','flat');hold on
for i=1:8
    b1.CData(i,:) = [55,140,230]/255;
end
b1.CData(9,:) = [237,145,29]/255;
ylabel('ppm','FontName','Arial','FontSize',14);
xlim([0.25,9.75])
ylim([15,20])
set(gca,'YTick', [15:1:20]);
set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2015:1:2023],'FontName','Arial','fontsize',14)
xlabel('Year','FontName','Arial','FontSize',14)
title('Difference in CO_2 mole fraction between April and August in East Trout Lake (ETL) site','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.057257909159624 1.10148392156519 0],'FontName','Arial','FontSize',14,'fontweight','bold')

nexttile([2,2])
yyaxis left
x=1:length(co2_2023);
p1=plot(co2_2023,'-','LineWidth',2,'MarkerSize',15,'color','r'); hold on
p2=plot(co2_mean,'-','LineWidth',2,'MarkerSize',15,'color','k');
p4=plot([0,10],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[157,157,157]/255);hold on
fill([x, fliplr(x)], [co2_mean, fliplr(co2_mean+co2_std)],'k','linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [co2_mean, fliplr(co2_mean-co2_std)],'k','linestyle', 'none', 'FaceAlpha',0.2);
ylim([-16,16])
set(gca,'YTick',[-16:8:16],'FontName','Arial','fontsize',14)
set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'M','A','M','J','J','A','S','O','N'},'FontName','Arial','fontsize',14)
xlabel('Month','FontName','Arial','FontSize',14)
ylabel('Detrended CO_2 mole fraction (ppm)','FontName','Arial','FontSize',14)
title('ETL','FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'YColor','k')
text('string','b','Units','normalized','position',[-0.146321509184895 1.10148392156519 0],'FontName','Arial','FontSize',14,'fontweight','bold')

yyaxis right
p3=plot(NEP_anomaly_mean,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255); hold on
fill([x, fliplr(x)], [NEP_anomaly_mean, fliplr(NEP_anomaly_mean+NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [NEP_anomaly_mean, fliplr(NEP_anomaly_mean-NEP_anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
ylim([-50,50])
set(gca,'YTick',[-50:25:50],'FontName','Arial','fontsize',14)
ylabel('NEP anomalies','FontName','Arial','FontSize',14)
xlim([0.5,9.5])
legend([p1 p2,p3],{'2023','Average level','NEP anomalies'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
set(gca,'YColor',[44,107,179]/255)


nexttile([2,2])
b=bar(co2_anomaly,'FaceColor',[153,221,227]/255); hold on
ylabel('ppm','FontName','Arial','FontSize',14);
xlabel('Season','FontName','Arial','FontSize',14)
ylim([-0.6,1.2])
set(gca,'YTick',[-0.6:0.6:1.2],'FontName','Arial','fontsize',14)

set(gca,'XTick',[1:1:3],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'MAM','JJA','SON'},'FontName','Arial','fontsize',14)
title('Detrended CO_2 mole fraction anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.146321509184895 1.10148392156519 0],'FontName','Arial','FontSize',14,'fontweight','bold')

set(gcf,'unit','centimeters','position',[26.431875000000005,11.906250000000002,28.17812500000001,20.82270833333334]);
result=['E:\phd_file\Boreal_North_America\Result\V6\site_CO2_line.png']
% print(result,f,'-r600','-dpng');
 
