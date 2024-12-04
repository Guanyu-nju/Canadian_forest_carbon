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
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
for year=2015:2023

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\year\NEE_' num2str(year) '.tif']);
    NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;
    NEE_list(year-2014)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));

    BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_annual\Opt_NEE_BEPS_GFAS_' num2str(year) '.tif']);
    BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_annual\Opt_NEE_BEPS_GFED_' num2str(year) '.tif']);
    CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_annual\Opt_NEE_CASA_GFAS_' num2str(year) '.tif']);
    CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_annual\Opt_NEE_CASA_GFED_' num2str(year) '.tif']);
    NEE_std_list(1,year-2014)=nansum(nansum(BEPS_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_std_list(2,year-2014)=nansum(nansum(BEPS_GFED_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_std_list(3,year-2014)=nansum(nansum(CASA_GFAS_NEE.*pixel_mask.*area_grid/(10^15)));
    NEE_std_list(4,year-2014)=nansum(nansum(CASA_GFED_NEE.*pixel_mask.*area_grid/(10^15)));

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GPP_list(year-2014)=nansum(nansum(GPP_mean_temp.*pixel_mask.*area_grid/(10^15)));

    GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\yearly\GOSIF_GPP_' num2str(year) '.tif']);
    Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\yearly\FluxSat_GPP_' num2str(year) '.tif']);

    GPP_std_list(1,year-2014)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
    GPP_std_list(2,year-2014)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\year\ER_' num2str(year) '.tif']);
    ER_list(year-2014)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));


end
NEE_mean=nanmean(NEE_mean,3);
NEE_std=nanstd(NEE_std_list);

for i=1:length(GPP_std_list)
    GPP_temp=GPP_std_list(:,i);
    GPP_std(i)=std(bootstrp(1000,@mean,GPP_temp));
end

ER_std=sqrt(power(NEE_std,2)+power(GPP_std,2));
    
total_list=[GPP_list',ER_list'];
total_std_list=[GPP_std',ER_std'];
%% calculate mean value (GCB2024)
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
for year=2015:2023

    GCB_NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\year\NEE_' num2str(year) '.tif']);
    GCB_NEE_mean(:,:,year-2014)=GCB_NNE_mean_temp.*pixel_mask;
    GCB_NEE_list(year-2014)=nansum(nansum(GCB_NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));

     A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']);


    GCB_NEE_std_list(1,year-2014)=nansum(nansum(A_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_std_list(2,year-2014)=nansum(nansum(B_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_std_list(3,year-2014)=nansum(nansum(C_NEE.*pixel_mask.*area_grid/(10^15)));
    GCB_NEE_std_list(4,year-2014)=nansum(nansum(D_NEE.*pixel_mask.*area_grid/(10^15)));
    
    % 
    % 
    % GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    % GPP_list(year-2014)=nansum(nansum(GPP_mean_temp.*pixel_mask.*area_grid/(10^15)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\year\ER_' num2str(year) '.tif']);
    GCB_ER_list(year-2014)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

   
end
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);
GCB_NEE_std=nanstd(GCB_NEE_std_list);
GCB_ER_std=sqrt(power(GCB_NEE_std,2)+power(GPP_std,2));

GCB_total_list=[GPP_list',GCB_ER_list'];
GCB_total_std_list=[GPP_std',GCB_ER_std'];
%%
% 面积占比
count_sum=sum(sum(~isnan(GCB_NEE_mean)));
count1=sum(sum(GCB_NEE_mean<0));
result=count1/(count_sum)
% NEE均值的标准差
nanstd(mean(NEE_std_list,2))
nanstd(mean(GCB_NEE_std_list,2))
% GPP均值的标准差
std(bootstrp(1000,@mean,mean(GPP_std_list,2)))
% sqrt(power(nanstd(mean(NEE_std_list,2)),2)+power(std(bootstrp(1000,@mean,mean(GPP_std_list,2))),2))
% sqrt(power(nanstd(mean(GCB_NEE_std_list,2)),2)+power(std(bootstrp(1000,@mean,mean(GPP_std_list,2))),2))
% ER均值的标准差
sqrt(power(0.08,2)+power(0.24,2))
sqrt(power(0.1,2)+power(0.24,2))
% 2023NEP变化标准差
std(NEE_std_list(:,end)-mean(NEE_std_list,2))
% 2023ER变化
ER_list(end)-mean(ER_list)
GCB_ER_list(end)-mean(GCB_ER_list)
% 2023GPP变化
GPP_list(end)-mean(GPP_list)
% 2023GPP变化标准差
std(bootstrp(1000,@mean,GPP_std_list(:,end)-mean(GPP_std_list,2)))
% 2023ER变化标准差
sqrt(power(std(NEE_std_list(:,end)-mean(NEE_std_list,2)),2)+power(std(bootstrp(1000,@mean,GPP_std_list(:,end)-mean(GPP_std_list,2))),2))
sqrt(power(std(GCB_NEE_std_list(:,end)-mean(GCB_NEE_std_list,2)),2)+power(std(bootstrp(1000,@mean,GPP_std_list(:,end)-mean(GPP_std_list,2))),2))

% GCAS和GCB2024整体结果
nanmean([NEE_std_list;GCB_NEE_std_list])
nanstd([NEE_std_list;GCB_NEE_std_list])
%%
f=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% part 1 *********************************************************************************
nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,-NEE_mean,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(nclCM(399));
caxis([-120,120]);
title('NEP (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string',['gC m^{-2} yr^{-1}'],'fontsize',14);
text('string','a','Units','normalized','position',[-0.0492086441485542 1.13824862744754 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 2 *********************************************************************************
nexttile
x=1:length(NEE_list);
yyaxis left
plot([0,20],[nanmean(-NEE_list),nanmean(-NEE_list)],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
p1=plot(-NEE_list,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
fill([x, fliplr(x)], [-NEE_list, fliplr(-NEE_list+NEE_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [-NEE_list, fliplr(-NEE_list-NEE_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
ylim([-0.3,0.6])
set(gca,'YTick', [-0.3:0.3:0.6]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,9.5])
set(gca,'YColor',[44,107,179]/255)

yyaxis right

p2=plot([0,20],[nanmean(total_list(:,1)),nanmean(total_list(:,1))],'--','LineWidth',1,'color',[17,119,51]/255);hold on
p3=plot([0,20],[nanmean(total_list(:,2)),nanmean(total_list(:,2))],'--','LineWidth',1,'color',[193,0,1]/255);

b=bar(total_list,'FaceColor','flat'); hold on
b(1).CData = [17,119,51]/255;
b(2).CData = [193,0,1]/255;

total_list(end,1)-nanmean(total_list(:,1))

% %分组误差棒
[M,N]=size(total_std_list);
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),total_list(:,:),total_std_list(:,:), ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

ylim([2.5,6.5])
set(gca,'YTick', [2:1:6]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);

set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2015:2023],'FontName','Arial','fontsize',14)
text('string','b','Units','normalized','position',[-0.135321692659772 1.04590643920493 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 b],{'NEP','GPP','ER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','off','Location','north')
xlabel('Year','FontName','Arial','FontSize',14)

% part 3 *********************************************************************************
nexttile
m_proj('lambert','long',[-150 -50],'lat',[42 75]);
m_patch(bou_canX,bou_canY,[225 225 227]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,-GCB_NEE_mean,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[48:10:78],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[48:10:78]); %对齐格网
colormap(nclCM(399));
caxis([-120,120]);
title('NEP (GCB2024)','FontName','Arial','FontSize',14,'fontweight','bold')
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string',['gC m^{-2} yr^{-1}'],'fontsize',14);
text('string','c','Units','normalized','position',[-0.0492086441485542 1.13824862744754 0],'FontName','Arial','FontSize',14,'fontweight','bold')

% part 4 *********************************************************************************
nexttile
x=1:length(NEE_list);
yyaxis left
plot([0,20],[nanmean(-GCB_NEE_list),nanmean(-GCB_NEE_list)],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
p1=plot(-GCB_NEE_list,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
fill([x, fliplr(x)], [-GCB_NEE_list, fliplr(-GCB_NEE_list+GCB_NEE_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [-GCB_NEE_list, fliplr(-GCB_NEE_list-GCB_NEE_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
ylim([-0.3,0.6])
set(gca,'YTick', [-0.3:0.3:0.6]);
set(gca,'YColor','k')
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,9.5])
set(gca,'YColor',[44,107,179]/255)

yyaxis right

p2=plot([0,20],[nanmean(GCB_total_list(:,1)),nanmean(GCB_total_list(:,1))],'--','LineWidth',1,'color',[17,119,51]/255);hold on
p3=plot([0,20],[nanmean(GCB_total_list(:,2)),nanmean(GCB_total_list(:,2))],'--','LineWidth',1,'color',[193,0,1]/255);

b=bar(GCB_total_list,'FaceColor','flat'); hold on
b(1).CData = [17,119,51]/255;
b(2).CData = [193,0,1]/255;
% GCB_total_list(end,1)-nanmean(GCB_total_list(:,1))

% %分组误差棒
[M,N]=size(total_std_list);
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),GCB_total_list(:,:),GCB_total_std_list(:,:), ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

ylim([2.5,6.5])
set(gca,'YTick', [2:1:6]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);

set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',[2015:2023],'FontName','Arial','fontsize',14)
text('string','d','Units','normalized','position',[-0.135321692659772 1.04590643920493 0],'FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 b],{'NEP','GPP','ER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','off','Location','north')
xlabel('Year','FontName','Arial','FontSize',14)


set(gcf,'unit','centimeters','position',[14.605000000000002,6.2865,34.798,20.277666666666672]);
result=['E:\phd_file\Boreal_North_America\Result\V6\NEP_map_line.png']
% print(result,f,'-r600','-dpng' );