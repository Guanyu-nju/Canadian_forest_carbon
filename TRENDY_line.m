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

        % 计算 NEP (逐月)
        annual_NEP_combined(:, :, :, y, f) = annual_GPP_combined(:, :, :, y, f) - annual_ER_combined(:, :, :, y, f);

        % 计算 ER (年度总量)
        annual_ER_sum(:, :, y, f) = annual_ra_sum(:, :, y, f) + annual_rh_sum(:, :, y, f);

        % 计算 NEP (年度总量)
        annual_NEP_sum(:, :, y, f) = annual_GPP_sum(:, :, y, f) - annual_ER_sum(:, :, y, f);
    end
end

% 逐年计算均值
for y = 1:size(annual_ER_sum,3)
    % 计算 NEP 的均值
    mean_NEP(:, :, y) = nanmean(annual_NEP_sum(:, :, y, :), 4);
    % 计算 ER 的均值
    mean_ER(:, :, y) = nanmean(annual_ER_sum(:, :, y, :), 4);
    % 计算 GPP 的均值
    mean_GPP(:, :, y) = nanmean(annual_GPP_sum(:, :, y, :), 4);
end
for y = 1:size(annual_ER_sum,3)
    % NEP 总量
    NEP_mean(1, y) = nansum(nansum(mean_NEP(:, :, y) .* area_grid))/(10^15);
    % ER 总量
    ER_mean(1, y) = nansum(nansum(mean_ER(:, :, y) .* area_grid))/(10^15);
    % GPP 总量
    GPP_mean(1, y) = nansum(nansum(mean_GPP(:, :, y) .* area_grid))/(10^15);
end

% 遍历每个nc文件结果
for f = 1:size(annual_ER_sum,4)
    for y = 1:size(annual_ER_sum,3)
        % 提取每年的二维矩阵
        NEP_yearly = annual_NEP_sum(:, :, y, f); % NEP
        ER_yearly = annual_ER_sum(:, :, y, f);   % ER
        GPP_yearly = annual_GPP_sum(:, :, y, f); % GPP

        % 加权求和：矩阵逐元素相乘并求和
        total_NEP_list(f, y) = nansum(nansum(NEP_yearly .* area_grid))/(10^15);
        total_ER_list(f, y) = nansum(nansum(ER_yearly .* area_grid))/(10^15);
        total_GPP_list(f, y) = nansum(nansum(GPP_yearly .* area_grid))/(10^15);
    end
end
total_NEP_std=nanstd(total_NEP_list);
total_ER_std=nanstd(total_ER_list);
total_GPP_std=nanstd(total_GPP_list);
%%
ff=figure
t = tiledlayout(1,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
x=1:length(NEP_mean);
yyaxis left
plot([0.5,9.5],[nanmean(NEP_mean),nanmean(NEP_mean)],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
p1=plot(NEP_mean,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
fill([x, fliplr(x)], [NEP_mean, fliplr(NEP_mean+total_NEP_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [NEP_mean, fliplr(NEP_mean-total_NEP_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
ylim([-0.1,0.5])
% set(gca,'YTick', [-0.3:0.3:0.6]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,9.5])
set(gca,'YColor',[44,107,179]/255)

yyaxis right

p2=plot([0.5,9.5],[nanmean(GPP_mean),nanmean(GPP_mean)],'--','LineWidth',1,'color',[17,119,51]/255);hold on
p3=plot([0.5,9.5],[nanmean(ER_mean),nanmean(ER_mean)],'--','LineWidth',1,'color',[193,0,1]/255);

b=bar([GPP_mean;ER_mean]','FaceColor','flat'); hold on
b(1).CData = [17,119,51]/255;
b(2).CData = [193,0,1]/255;


% %分组误差棒
[M,N]=size([GPP_mean;ER_mean]');
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),[GPP_mean;ER_mean]',[total_GPP_std;total_ER_std]', ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

ylim([3.7,8])
set(gca,'YTick', [4:1:8]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);

set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2015:2023],'FontName','Arial','fontsize',14)
title('Ensemble mean of TRENDY DGVMs','FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 b],{'NEP','GPP','ER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','off','Location','north')
xlabel('Year','FontName','Arial','FontSize',14)
set(gcf,'unit','centimeters','position',[26.431875000000005,21.642916666666668,16.80104166666667,11.08604166666667]);
result=['E:\phd_file\Boreal_North_America\Result\V6\TRENDY_NEP_line.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCN","ORCHIDEE","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';


for i=1:length(filelist)

    nexttile
    x=1:length(total_NEP_list(i,:));
    yyaxis left
    plot([0.5,9.5],[nanmean(total_NEP_list(i,:)),nanmean(total_NEP_list(i,:))],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
    p1=plot(total_NEP_list(i,:),'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
    if i==1
        ylim([-0.05,0.3])
        set(gca,'YTick', [0:0.1:0.3]);
    elseif i==2
        ylim([-0.15,0.2])
        set(gca,'YTick', [-0.1:0.1:0.2]);
    elseif i==3
        ylim([0.05,0.5])
        set(gca,'YTick', [0.1:0.1:0.5]);
    elseif i==4
        ylim([0.05,0.17])
        set(gca,'YTick', [0.05:0.03:0.17]);
    elseif i==5
        ylim([0,0.4])
        set(gca,'YTick', [0.1:0.1:0.5]);
    elseif i==6
        ylim([-0.1,0.25])
        set(gca,'YTick', [-0.1:0.1:0.2]);
    elseif i==7
        ylim([0.1,0.5])
        set(gca,'YTick', [0.1:0.1:0.5]);
    elseif i==8
        ylim([0.15,0.4])
        set(gca,'YTick', [0.1:0.1:0.4]);
    elseif i==9
        ylim([-0.05,0.42])
        set(gca,'YTick', [0:0.1:0.4]);
    elseif i==10
        ylim([-0.1,0.4])
        set(gca,'YTick', [-0.1:0.1:0.4]);
    elseif i==11
        ylim([-0.05,0.4])
        set(gca,'YTick', [-0.1:0.1:0.4]);
    elseif i==12
        ylim([0.05,0.3])
        set(gca,'YTick', [0.1:0.1:0.3]);
    elseif i==13
        ylim([0.05,0.4])
        set(gca,'YTick', [0.1:0.1:0.4]);
    elseif i==14
        ylim([0.1,0.55])
        set(gca,'YTick', [0.1:0.1:0.5]);
    end

    ylabel('PgC yr^{-1}','FontName','Arial','FontSize',12);
    xlim([0.5,9.5])
    set(gca,'YColor',[44,107,179]/255)

    yyaxis right

    p2=plot([0.5,9.5],[nanmean(total_GPP_list(i,:)),nanmean(total_GPP_list(i,:))],'--','LineWidth',1,'color',[17,119,51]/255);hold on
    p3=plot([0.5,9.5],[nanmean(total_ER_list(i,:)),nanmean(total_ER_list(i,:))],'--','LineWidth',1,'color',[193,0,1]/255);

    b=bar([total_GPP_list(i,:);total_ER_list(i,:)]','FaceColor','flat'); hold on
    b(1).CData = [17,119,51]/255;
    b(2).CData = [193,0,1]/255;
    if i==1
        ylim([3,7.4])
        set(gca,'YTick', [3:1:8]);
    elseif i==2
        ylim([2.8,8])
        set(gca,'YTick', [3:1:8]);
    elseif i==3
        ylim([3.8,11.5])
        set(gca,'YTick', [3:1:11]);
    elseif i==4
        ylim([2.5,6.5])
        set(gca,'YTick', [3:1:6]);
    elseif i==5
        ylim([3.3,11])
        set(gca,'YTick', [3:1:11]);
    elseif i==6
        ylim([3.8,11.5])
        set(gca,'YTick', [3:1:11]);
    elseif i==7
        ylim([1.5,5])
        set(gca,'YTick', [1:1:5]);
    elseif i==8
        ylim([1.5,12])
        set(gca,'YTick', [1:1:9]);
    elseif i==9
        ylim([4.5,11.5])
        set(gca,'YTick', [3:1:11]);
    elseif i==10
        ylim([2.5,7.4])
        set(gca,'YTick', [2:1:7]);
    elseif i==11
        ylim([2.5,10])
        set(gca,'YTick', [3:1:11]);
    elseif i==12
        ylim([2.6,10])
        set(gca,'YTick', [3:1:11]);
    elseif i==13
        ylim([2.2,6])
        set(gca,'YTick', [3:1:11]);
    elseif i==14
        ylim([3,11])
        set(gca,'YTick', [3:1:11]);

    end
    ylabel('PgC yr^{-1}','FontName','Arial','FontSize',12);

    set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',12)
    set(gca,'XTickLabel',[2015:2023],'FontName','Arial','fontsize',12)
    text('string',panellist{i},'Units','normalized','position',[-0.152748545258527 1.06288467740369 0],'FontName','Arial','FontSize',12,'fontweight','bold')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')

end
lgd = legend([p1 b],{'NEP','GPP','ER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','off');
% lgd.Layout.Tile = 15;
set(gcf,'unit','centimeters','position',[4.048125000000001,2.037291666666667,57.99666666666668,32.385000000000005]);
result=['E:\phd_file\Boreal_North_America\Result\V6\TRENDY_indivual_NEP_line.png']
% print(result,ff,'-r600','-dpng' );

