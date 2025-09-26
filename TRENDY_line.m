%% Mask create
clc
clear
% Canada mask
pixel_mask=importdata('E:\phd_file\Boreal_North_America\Canada_mask.tif');
forest_mask=importdata("E:\phd_file\Boreal_North_America\region_lu.tif");
forest_mask(forest_mask~=1)=nan;
pixel_mask=pixel_mask.*forest_mask;

% Create environment
data = geotiffread('E:\phd_file\Boreal_North_America\region_lu.tif');
info = geotiffinfo('E:\phd_file\Boreal_North_America\region_lu.tif');
[m, n] = size(data);

% Extract latitude coordinates from the first column
k = 1;
for i = 1:m
    for j = 1:1
        [lat, lon] = pix2latlon(info.RefMatrix, i, j);   % Read latitude for all rows in first column
        lat_(k, :) = lat; % Store latitude data as a column
        k = k + 1;
    end
end

% Extract longitude coordinates from the first row
k = 1;
for ii = 1:1
    for jj = 1:n
        [lat, lon] = pix2latlon(info.RefMatrix, ii, jj);   % Read longitude for all columns in first row
        lon_(k, :) = lon;  % Store longitude data as a column
        k = k + 1;
    end
end

% Create coordinate grids using meshgrid
[lon1, lat1] = meshgrid(lon_, lat_);

Boundry = shaperead("E:\phd_file\yuling_shiliang\countries.shp");
bou_canX = [Boundry(:).X];
bou_canY= [Boundry(:).Y];

clearvars -except lon1 lat1 pixel_mask bou_canX bou_canY m n
%%
load TRENDY_variable.mat

area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;


% calculate ER NEP
for f = 1:size(annual_GPP_combined,5)
    for y = 1:size(annual_GPP_combined,4)
        % calculate ER (monthly)
        annual_ER_combined(:, :, :, y, f) = annual_ra_combined(:, :, :, y, f) + annual_rh_combined(:, :, :, y, f);

        % calculate NEP (monthly)
        annual_NEP_combined(:, :, :, y, f) = annual_GPP_combined(:, :, :, y, f) - annual_ER_combined(:, :, :, y, f);

        % calculate ER (annual)
        annual_ER_sum(:, :, y, f) = annual_ra_sum(:, :, y, f) + annual_rh_sum(:, :, y, f);

        % calculate NEP (annual)
        annual_NEP_sum(:, :, y, f) = annual_GPP_sum(:, :, y, f) - annual_ER_sum(:, :, y, f);
    end
end

% calculate mean value
for y = 1:size(annual_ER_sum,3)
    % /NEP
    mean_NEP(:, :, y) = nanmean(annual_NEP_sum(:, :, y, :), 4);
    % ER
    mean_ER(:, :, y) = nanmean(annual_ER_sum(:, :, y, :), 4);
    %  GPP
    mean_GPP(:, :, y) = nanmean(annual_GPP_sum(:, :, y, :), 4);
    % RZSM
    mean_RZSM(:, :, y) = nanmean(annual_mrso_sum(:, :, y, :), 4);
    % TEM
    mean_TEM(:, :, y) = nanmean(annual_tas_sum(:, :, y, :), 4);

end

% 
mean_NEP_2023anomalies=mean_NEP(:,:,end)-nanmean(mean_NEP(:,:,1:end-1),3);
mean_ER_2023anomalies=mean_ER(:,:,end)-nanmean(mean_ER(:,:,1:end-1),3);
mean_GPP_2023anomalies=mean_GPP(:,:,end)-nanmean(mean_GPP(:,:,1:end-1),3);

% seasonal anomaly

for f = 1:size(annual_ER_combined,5)
    for y = 1:size(annual_ER_combined,4)
       
        ER_year = annual_ER_combined(:, :, :, y, f); %  6-8 month
        GPP_year = annual_GPP_combined(:, :, :, y, f);
        NEP_year = annual_NEP_combined(:, :, :, y, f);

        for month=1:12
            ER_month_list(y,month,f)=nansum(nansum(ER_year(:,:,month).*area_grid/(10^15)));
            NEP_month_list(y,month,f)=nansum(nansum(NEP_year(:,:,month).*area_grid/(10^15)));
            GPP_month_list(y,month,f)=nansum(nansum(GPP_year(:,:,month).*area_grid/(10^15)));
        end

    end
end
% 
for f=1:size(ER_month_list,3)
    ER_month_temp=ER_month_list(:,:,f);
    ER_month_2023anomaly(f,:)=ER_month_temp(end,:)-nanmean(ER_month_temp(1:end-1,:));

    NEP_month_temp=NEP_month_list(:,:,f);
    NEP_month_2023anomaly(f,:)=NEP_month_temp(end,:)-nanmean(NEP_month_temp(1:end-1,:));

    GPP_month_temp=GPP_month_list(:,:,f);
    GPP_month_2023anomaly(f,:)=GPP_month_temp(end,:)-nanmean(GPP_month_temp(1:end-1,:));
end

mean_ER_month_2023anomaly=nanmean(ER_month_2023anomaly);
mean_ER_month_2023anomaly_std=nanstd(ER_month_2023anomaly);

mean_NEP_month_2023anomaly=nanmean(NEP_month_2023anomaly);
mean_NEP_month_2023anomaly_std=nanstd(NEP_month_2023anomaly);

mean_GPP_month_2023anomaly=nanmean(GPP_month_2023anomaly);
mean_GPP_month_2023anomaly_std=nanstd(GPP_month_2023anomaly);


for y = 1:size(annual_ER_sum,3)
    % NEP
    NEP_mean(1, y) = nansum(nansum(mean_NEP(:, :, y) .* area_grid))/(10^15);
    % ER
    ER_mean(1, y) = nansum(nansum(mean_ER(:, :, y) .* area_grid))/(10^15);
    % GPP
    GPP_mean(1, y) = nansum(nansum(mean_GPP(:, :, y) .* area_grid))/(10^15);

    % RZSM
    RZSM_mean(1, y) = nansum(nansum(mean_RZSM(:, :, y).*area_grid))/(nansum(nansum(area_grid)));
    % TEM
    TEM_mean(1, y) = nansum(nansum(mean_TEM(:, :, y).*area_grid))/(nansum(nansum(area_grid)));
end

% 
for f = 1:size(annual_ER_sum,4)
    for y = 1:size(annual_ER_sum,3)
        % 
        NEP_yearly = annual_NEP_sum(:, :, y, f); % NEP
        ER_yearly = annual_ER_sum(:, :, y, f);   % ER
        GPP_yearly = annual_GPP_sum(:, :, y, f); % GPP
        TEM_yearly = annual_tas_sum(:, :, y, f); % TEM
        RZSM_yearly = annual_mrso_sum(:, :, y, f); % RZSM

        % 
        total_NEP_list(f, y) = nansum(nansum(NEP_yearly .* area_grid))/(10^15);
        total_ER_list(f, y) = nansum(nansum(ER_yearly .* area_grid))/(10^15);
        total_GPP_list(f, y) = nansum(nansum(GPP_yearly .* area_grid))/(10^15);
        total_TEM_list(f, y) = nansum(nansum(TEM_yearly.*area_grid))/(nansum(nansum(area_grid)));
        total_RZSM_list(f, y) = nansum(nansum(RZSM_yearly.*area_grid))/(nansum(nansum(area_grid)));

    end
end
total_NEP_std=nanstd(total_NEP_list);
total_ER_std=nanstd(total_ER_list);
total_GPP_std=nanstd(total_GPP_list);

%%
ff=figure
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','centimeters','position',[16.774583333333336,6.773333333333334,35.956875,19.711458333333336]);


n1=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,mean_NEP_2023anomalies.*pixel_mask,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
c = redblue();
colormap(n1,nclCM(399));

caxis([-40,40]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.0917680225652848 1.1119605181233 0],'FontName','Arial','FontSize',20,'fontweight','bold')

n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,mean_GPP_2023anomalies.*pixel_mask,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('GPP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.0917680225652848 1.1119605181233 0],'FontName','Arial','FontSize',20,'fontweight','bold')

n3=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,mean_ER_2023anomalies.*pixel_mask,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
c = redblue();
colormap(n3,flipud(nclCM(399)));

caxis([-160,160]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('TER anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.0917680225652848 1.1119605181233 0],'FontName','Arial','FontSize',20,'fontweight','bold')

n4=nexttile
x=1:length(mean_NEP_month_2023anomaly);
plot([0,20],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[157,157,157]/255);hold on
p1=plot(mean_NEP_month_2023anomaly,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);
fill([x, fliplr(x)], [mean_NEP_month_2023anomaly, fliplr(mean_NEP_month_2023anomaly+mean_NEP_month_2023anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [mean_NEP_month_2023anomaly, fliplr(mean_NEP_month_2023anomaly-mean_NEP_month_2023anomaly_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p2=plot(mean_GPP_month_2023anomaly,'-','LineWidth',2,'MarkerSize',15,'color',[17,119,51]/255);
fill([x, fliplr(x)], [mean_GPP_month_2023anomaly, fliplr(mean_GPP_month_2023anomaly+mean_GPP_month_2023anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [mean_GPP_month_2023anomaly, fliplr(mean_GPP_month_2023anomaly-mean_GPP_month_2023anomaly_std)],[17,119,51]/255,'linestyle', 'none', 'FaceAlpha',0.2);
p3=plot(mean_ER_month_2023anomaly,'-','LineWidth',2,'MarkerSize',15,'color',[193,0,1]/255);hold on
fill([x, fliplr(x)], [mean_ER_month_2023anomaly, fliplr(mean_ER_month_2023anomaly+mean_ER_month_2023anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [mean_ER_month_2023anomaly, fliplr(mean_ER_month_2023anomaly-mean_ER_month_2023anomaly_std)],[193,0,1]/255,'linestyle', 'none', 'FaceAlpha',0.2);


ylim([-0.05,0.17])
% set(gca,'YTick', [-0.1:0.1:0.2]);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
text('string','d','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
% title('Monthly flux anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_NEP_map_line.png']
% print(result,ff,'-r600','-dpng' );
%% NEP map

filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[135,193,2178,1024]);

for y = 1:size(annual_ER_sum,4)
    
    annual_NEP_temp=annual_NEP_sum(:, :, :, y);
    annual_NEP_2023_anomaly_temp=annual_NEP_temp(:,:,end)-nanmean(annual_NEP_temp(:,:,1:end-1),3);

    nexttile
    m_proj('lambert','long',[-150 -50],'lat',[41 75]);
    m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
    % m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
    m_pcolor(lon1,lat1,annual_NEP_2023_anomaly_temp.*pixel_mask,'linestyle','none');
    % m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
    m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
        'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
    colormap(nclCM(399));
    caxis([-60,60]);
    h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-60:30:60]);
    set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
    ss=filelist{y};
    ss=['NEP anomalies (', ss{1}, ')'];
    title(ss,'FontName','Arial','FontSize',14,'fontweight','bold')
    text('string',panellist{y},'Units','normalized','position',[-0.152748545258527 1.06288467740369 0],'FontName','Arial','FontSize',20,'fontweight','bold')

end
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_NEP_map.png']
% print(result,ff,'-r600','-dpng' );
%% GPP map

filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[135,193,2178,1024]);

for y = 1:size(annual_ER_sum,4)
    
    annual_GPP_temp=annual_GPP_sum(:, :, :, y);
    annual_GPP_2023_anomaly_temp=annual_GPP_temp(:,:,end)-nanmean(annual_GPP_temp(:,:,1:end-1),3);

    nexttile
    m_proj('lambert','long',[-150 -50],'lat',[41 75]);
    m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
    % m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
    m_pcolor(lon1,lat1,annual_GPP_2023_anomaly_temp.*pixel_mask,'linestyle','none');
    % m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
    m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
        'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
    colormap(nclCM(399));
    caxis([-160,160]);
    h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
    set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
    ss=filelist{y};
    ss=['GPP anomalies (', ss{1}, ')'];
    title(ss,'FontName','Arial','FontSize',14,'fontweight','bold')
    text('string',panellist{y},'Units','normalized','position',[-0.152748545258527 1.06288467740369 0],'FontName','Arial','FontSize',20,'fontweight','bold')

end
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_GPP_map.png']
% print(result,ff,'-r600','-dpng' );
%% TER map

filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[135,193,2178,1024]);

for y = 1:size(annual_ER_sum,4)
    
    annual_ER_temp=annual_ER_sum(:, :, :, y);
    annual_ER_2023_anomaly_temp=annual_ER_temp(:,:,end)-nanmean(annual_ER_temp(:,:,1:end-1),3);

    nexttile
    m_proj('lambert','long',[-150 -50],'lat',[41 75]);
    m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
    % m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
    m_pcolor(lon1,lat1,annual_ER_2023_anomaly_temp.*pixel_mask,'linestyle','none');
    % m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
    m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
        'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
    colormap(flipud(nclCM(399)));
    caxis([-160,160]);
    h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-160:80:160]);
    set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
    ss=filelist{y};
    ss=['TER anomalies (', ss{1}, ')'];
    title(ss,'FontName','Arial','FontSize',14,'fontweight','bold')
    text('string',panellist{y},'Units','normalized','position',[-0.152748545258527 1.06288467740369 0],'FontName','Arial','FontSize',20,'fontweight','bold')

end
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_TER_map.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','pixels','position',[310,83,1750,1215]);

for i=1:length(filelist)
    nexttile
    % x=1:length(mean_NEP_month_2023anomaly);
    plot([0,20],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[157,157,157]/255);hold on
    p1=plot(NEP_month_2023anomaly(i,:),'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);
    p2=plot(GPP_month_2023anomaly(i,:),'-','LineWidth',2,'MarkerSize',15,'color',[17,119,51]/255);
    p3=plot(ER_month_2023anomaly(i,:),'-','LineWidth',2,'MarkerSize',15,'color',[193,0,1]/255);hold on

    ylim([-0.08,0.24])
    set(gca,'YTick', [-0.08:0.08:0.24]);
    xlim([0.5,12.5])
    set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
    set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)

    title(filelist{i},'FontName','Arial','FontSize',14,'fontweight','bold')
    ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
    text('string',panellist{i},'Units','normalized','position',[-0.203953364535635 1.10939630635028 0],'FontName','Arial','FontSize',20,'fontweight','bold')
    xlabel('Month','FontName','Arial','FontSize',14)

end
% legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
lgd = legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','on');
lgd.Layout.Tile = 15;
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_line2.png']
% print(result,ff,'-r600','-dpng' );
%%
ff=figure
t = tiledlayout(1,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','centimeters','position',[16.774583333333336,13.096875000000002,19.711458333333336,13.387916666666671]);

x=1:length(NEP_mean);
yyaxis left
plot([0.5,9.5],[nanmean(NEP_mean(1:end-1)),nanmean(NEP_mean(1:end-1))],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
p1=plot(NEP_mean,'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
fill([x, fliplr(x)], [NEP_mean, fliplr(NEP_mean+total_NEP_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
fill([x, fliplr(x)], [NEP_mean, fliplr(NEP_mean-total_NEP_std)],[44,107,179]/255,'linestyle', 'none', 'FaceAlpha',0.2);
ylim([-0.1,0.5])
% set(gca,'YTick', [-0.3:0.3:0.6]);
ylabel('PgC yr^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,9.5])
set(gca,'YColor',[44,107,179]/255)

yyaxis right

p2=plot([0.5,9.5],[nanmean(GPP_mean(1:end-1)),nanmean(GPP_mean(1:end-1))],'--','LineWidth',1,'color',[17,119,51]/255);hold on
p3=plot([0.5,9.5],[nanmean(ER_mean(1:end-1)),nanmean(ER_mean(1:end-1))],'--','LineWidth',1,'color',[193,0,1]/255);

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
% title('Ensemble mean of TRENDY DGVMs','FontName','Arial','FontSize',14,'fontweight','bold')
legend([p1 b],{'NEP','GPP','TER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','off','Location','north')
xlabel('Year','FontName','Arial','FontSize',14)
% text('string','d','Units','normalized','position',[-0.131105496685367 1.02314472887019 0],'FontName','Arial','FontSize',20,'fontweight','bold')

result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_NEP_line.png']
% print(result,ff,'-r600','-dpng' );
%%
filelist={"CABLE-POP","CLASSIC","CLM6.0","EDv3","ELM","IBIS","ISBA-CTRIP","JSBACH","LPJmL","LPX-Bern","OCNv2","ORCHIDEEv3","SDGVM","VISIT"};
panellist={"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"};
ff=figure
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';


for i=1:length(filelist)

    nexttile
    x=1:length(total_NEP_list(i,:));
    yyaxis left
    plot([0.5,9.5],[nanmean(total_NEP_list(i,1:end-1)),nanmean(total_NEP_list(i,1:end-1))],'--','LineWidth',1,'MarkerSize',3,'color',[44,107,179]/255);hold on
    p1=plot(total_NEP_list(i,:),'-','LineWidth',2,'MarkerSize',15,'color',[44,107,179]/255);hold on
    if i==1
        ylim([-0.05,0.3])
        set(gca,'YTick', [0:0.1:0.3]);
    elseif i==2
        ylim([-0.1,0.2])
        set(gca,'YTick', [-0.1:0.1:0.2]);
    elseif i==3
        ylim([0,0.4])
        set(gca,'YTick', [0:0.1:0.4]);
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
        ylim([5,11.5])
        set(gca,'YTick', [5:2:11]);
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
    text('string',panellist{i},'Units','normalized','position',[-0.152748545258527 1.06288467740369 0],'FontName','Arial','FontSize',20,'fontweight','bold')
    title(filelist{i},'FontName','Arial','FontSize',12,'fontweight','bold')

end

    lgd = legend([p1 b],{'NEP','GPP','TER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','on');
lgd.Layout.Tile = 15;

set(gcf,'unit','centimeters','position',[4.048125000000001,2.037291666666667,57.99666666666668,32.385000000000005]);
result=['E:\phd_file\Boreal_North_America\Result\V9\TRENDY_indivual_NEP_line.png']
% print(result,ff,'-r600','-dpng' );
%%
for year=2015:2023

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\year\NEE_' num2str(year) '.tif']);
    GCAS_NEE_mean(year-2014)=nansum(nansum(NNE_mean_temp.*area_grid/(10^15)));

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GCAS_GPP_mean(year-2014)=nansum(nansum(GPP_mean_temp.*area_grid/(10^15)));

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\year\ER_' num2str(year) '.tif']);
    GCAS_ER_mean(year-2014)=nansum(nansum(ER_mean_temp.*area_grid/(10^15)));

end

for i=1:14

    [rho,pval] = corrcoef(-GCAS_NEE_mean,total_NEP_list(i,:));
    NEP_corr(1,i)=rho(1,2);
    NEP_pvalue(1,i)=pval(1,2);

    [rho,pval] = corrcoef(GCAS_ER_mean,total_ER_list(i,:));
    ER_corr(1,i)=rho(1,2);
    ER_pvalue(1,i)=pval(1,2);

    [rho,pval] = corrcoef(GCAS_GPP_mean,total_GPP_list(i,:));
    GPP_corr(1,i)=rho(1,2);
    GPP_pvalue(1,i)=pval(1,2);
end
[rho,pval] = corrcoef(-GCAS_NEE_mean,NEP_mean)
[rho,pval] = corrcoef(GCAS_ER_mean,ER_mean)
[rho,pval] = corrcoef(GCAS_GPP_mean,GPP_mean)

