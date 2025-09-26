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
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;
%% calculate landflux value (GCAS)

for year=2015:2023

    % original fire
    GFAS_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\total_carbon\year\Fire_' num2str(year) '.tif']);
    GFED_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\total_carbon\year\Fire_' num2str(year) '.tif']);
    Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\year\Fire_' num2str(year) '.tif']);
    % case NEE
    BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_total_annual\Opt_NEE_BEPS_GFAS_' num2str(year) '.tif']);
    BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_total_annual\Opt_NEE_BEPS_GFED_' num2str(year) '.tif']);
    CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_total_annual\Opt_NEE_CASA_GFAS_' num2str(year) '.tif']);
    CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_total_annual\Opt_NEE_CASA_GFED_' num2str(year) '.tif']);
    % original landflux
    BEPS_GFAS_landflux=BEPS_GFAS_NEE+GFAS_fire;
    BEPS_GFED_landflux=BEPS_GFED_NEE+GFED_fire;
    CASA_GFAS_landflux=CASA_GFAS_NEE+GFAS_fire;
    CASA_GFED_landflux=CASA_GFED_NEE+GFED_fire;
    % other case NEE
    BEPS_GFAS_towards_GFED_NEE=BEPS_GFAS_landflux-GFED_fire;
    BEPS_GFAS_towards_mean_NEE=BEPS_GFAS_landflux-Mean_fire;
    BEPS_GFED_towards_GFAS_NEE=BEPS_GFED_landflux-GFAS_fire;
    BEPS_GFED_towards_mean_NEE=BEPS_GFED_landflux-Mean_fire;
    CASA_GFAS_towards_GFED_NEE=CASA_GFAS_landflux-GFED_fire;
    CASA_GFAS_towards_mean_NEE=CASA_GFAS_landflux-Mean_fire;
    CASA_GFED_towards_GFAS_NEE=CASA_GFED_landflux-GFAS_fire;
    CASA_GFED_towards_mean_NEE=CASA_GFED_landflux-Mean_fire;
    % GPP
    GOSIF_GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\yearly\GOSIF_GPP_' num2str(year) '.tif']);
    FluxSat_GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\yearly\FluxSat_GPP_' num2str(year) '.tif']);
    Mean_GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    % Fire list
    GFAS_fire_list(year-2014)=nansum(nansum(GFAS_fire.*area_grid/(10^15)));
    GFED_fire_list(year-2014)=nansum(nansum(GFED_fire.*area_grid/(10^15)));
    Mean_fire_list(year-2014)=nansum(nansum(Mean_fire.*area_grid/(10^15)));
    % landflux list
    GCAS_landflux_list(1,year-2014)=nansum(nansum(BEPS_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(2,year-2014)=nansum(nansum(BEPS_GFED_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(3,year-2014)=nansum(nansum(CASA_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_list(4,year-2014)=nansum(nansum(CASA_GFED_landflux.*area_grid/(10^15)));
    
    % GFAS-NEE list
    GCAS_GFAS_NEE_list(1,year-2014)=nansum(nansum(BEPS_GFAS_NEE.*area_grid/(10^15)));
    GCAS_GFAS_NEE_list(2,year-2014)=nansum(nansum(CASA_GFAS_NEE.*area_grid/(10^15)));
    GCAS_GFAS_NEE_list(3,year-2014)=nansum(nansum(BEPS_GFED_towards_GFAS_NEE.*area_grid/(10^15)));
    GCAS_GFAS_NEE_list(4,year-2014)=nansum(nansum(CASA_GFED_towards_GFAS_NEE.*area_grid/(10^15)));
    % GFED-NEE list
    GCAS_GFED_NEE_list(1,year-2014)=nansum(nansum(BEPS_GFED_NEE.*area_grid/(10^15)));
    GCAS_GFED_NEE_list(2,year-2014)=nansum(nansum(CASA_GFED_NEE.*area_grid/(10^15)));
    GCAS_GFED_NEE_list(3,year-2014)=nansum(nansum(BEPS_GFAS_towards_GFED_NEE.*area_grid/(10^15)));
    GCAS_GFED_NEE_list(4,year-2014)=nansum(nansum(CASA_GFAS_towards_GFED_NEE.*area_grid/(10^15)));
    % Meanfire-NEE list
    GCAS_Mean_fire_NEE_list(1,year-2014)=nansum(nansum(BEPS_GFAS_towards_mean_NEE.*area_grid/(10^15)));
    GCAS_Mean_fire_NEE_list(2,year-2014)=nansum(nansum(BEPS_GFED_towards_mean_NEE.*area_grid/(10^15)));
    GCAS_Mean_fire_NEE_list(3,year-2014)=nansum(nansum(CASA_GFAS_towards_mean_NEE.*area_grid/(10^15)));
    GCAS_Mean_fire_NEE_list(4,year-2014)=nansum(nansum(CASA_GFED_towards_mean_NEE.*area_grid/(10^15)));
    % GPP list
    GOSIF_GPP_list(year-2014)=nansum(nansum(GOSIF_GPP.*area_grid/(10^15)));
    FluxSat_GPP_list(year-2014)=nansum(nansum(FluxSat_GPP.*area_grid/(10^15)));
    Mean_GPP_list(year-2014)=nansum(nansum(Mean_GPP.*area_grid/(10^15)));

end

% ER list
GCAS_GFAS_GOSIF_ER_list=mean(GCAS_GFAS_NEE_list)+GOSIF_GPP_list;
GCAS_GFAS_FluxSat_ER_list=mean(GCAS_GFAS_NEE_list)+FluxSat_GPP_list;
GCAS_GFAS_Mean_GPP_ER_list=mean(GCAS_GFAS_NEE_list)+Mean_GPP_list;
GCAS_GFED_GOSIF_ER_list=mean(GCAS_GFED_NEE_list)+GOSIF_GPP_list;
GCAS_GFED_FluxSat_ER_list=mean(GCAS_GFED_NEE_list)+FluxSat_GPP_list;
GCAS_GFED_Mean_GPP_ER_list=mean(GCAS_GFED_NEE_list)+Mean_GPP_list;
GCAS_Mean_fire_GOSIF_ER_list=mean(GCAS_Mean_fire_NEE_list)+GOSIF_GPP_list;
GCAS_Mean_fire_FluxSat_ER_list=mean(GCAS_Mean_fire_NEE_list)+FluxSat_GPP_list;
GCAS_Mean_fire_Mean_GPP_ER_list=mean(GCAS_Mean_fire_NEE_list)+Mean_GPP_list;


% 2023 Flux value
GCAS_landflux_2023=mean(GCAS_landflux_list(:,end));
GCAS_Byrne_NEE_2023=GCAS_landflux_2023-0.647;
GCAS_GFAS_NEE_2023mean=mean(GCAS_GFAS_NEE_list(:,end));
GCAS_GFED_NEE_2023mean=mean(GCAS_GFED_NEE_list(:,end));
GCAS_Mean_fire_NEE_2023mean=mean(GCAS_Mean_fire_NEE_list(:,end));
Mean_GPP_2023=Mean_GPP_list(end);
GCAS_ER_2023=GCAS_Mean_fire_Mean_GPP_ER_list(end);
Mean_fire_2023=Mean_fire_list(end);

GCAS_GFAS_GOSIF_ER_2023=GCAS_GFAS_GOSIF_ER_list(end);
GCAS_GFAS_FluxSat_ER_2023=GCAS_GFAS_FluxSat_ER_list(end);
GCAS_GFAS_Mean_GPP_ER_2023=GCAS_GFAS_Mean_GPP_ER_list(end);
GCAS_GFED_GOSIF_ER_2023=GCAS_GFED_GOSIF_ER_list(end);
GCAS_GFED_FluxSat_ER_2023=GCAS_GFED_FluxSat_ER_list(end);
GCAS_GFED_Mean_GPP_ER_2023=GCAS_GFED_Mean_GPP_ER_list(end);
GCAS_Mean_fire_GOSIF_ER_2023=GCAS_Mean_fire_GOSIF_ER_list(end);
GCAS_Mean_fire_FluxSat_ER_2023=GCAS_Mean_fire_FluxSat_ER_list(end);
GCAS_Bynre_GOSIF_ER_2023=GCAS_Byrne_NEE_2023+GOSIF_GPP_list(end);
GCAS_Bynre_FluxSat_ER_2023=GCAS_Byrne_NEE_2023+FluxSat_GPP_list(end);
GCAS_Bynre_Mean_GPP_ER_2023=GCAS_Byrne_NEE_2023+Mean_GPP_2023;



GCAS_landflux_2023std=std(GCAS_landflux_list(:,end));
GCAS_GFAS_NEE_2023std=sqrt(power(GFAS_fire_list(end)*0.2,2)+power(GCAS_landflux_2023std,2));
GCAS_GFED_NEE_2023std=sqrt(power(GFED_fire_list(end)*0.2,2)+power(GCAS_landflux_2023std,2));
GCAS_Bynre_NEE_2023std=sqrt(power(0.647*0.2,2)+power(GCAS_landflux_2023std,2));
GCAS_Mean_fire_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(GCAS_landflux_2023std,2));
GPP_2023std=std(bootstrp(1000,@mean,[GOSIF_GPP_list(end),FluxSat_GPP_list(end)]));
Fire_2023std=Mean_fire_list(end)*0.2;
GCAS_ER_2023std=sqrt(power(Fire_2023std,2)+power(GCAS_landflux_2023std,2)+power(GPP_2023std,2));

% 2023 Anomaly relative to 2015-2022
GCAS_landflux_2023anomaly_mean=mean(GCAS_landflux_list(:,end))-mean(mean(GCAS_landflux_list(:,1:end-1)));
GCAS_GFAS_NEE_2023anomaly_mean=mean(GCAS_GFAS_NEE_list(:,end))-mean(mean(GCAS_GFAS_NEE_list(:,1:end-1)));
GCAS_GFED_NEE_2023anomaly_mean=mean(GCAS_GFED_NEE_list(:,end))-mean(mean(GCAS_GFED_NEE_list(:,1:end-1)));
GCAS_Mean_fire_NEE_2023anomaly_mean=mean(GCAS_Mean_fire_NEE_list(:,end))-mean(mean(GCAS_Mean_fire_NEE_list(:,1:end-1)));
Mean_GPP_2023anomaly=Mean_GPP_list(end)-mean(Mean_GPP_list(1:end-1));
GCAS_ER_2023anomaly=GCAS_Mean_fire_Mean_GPP_ER_list(end)-mean(GCAS_Mean_fire_Mean_GPP_ER_list(1:end-1));
Mean_fire_2023anomaly=Mean_fire_list(end)-mean(Mean_fire_list(1:end-1));

GCAS_GFAS_GOSIF_ER_2023anomaly=GCAS_GFAS_GOSIF_ER_list(end)-mean(GCAS_GFAS_GOSIF_ER_list(1:end-1));
GCAS_GFAS_FluxSat_ER_2023anomaly=GCAS_GFAS_FluxSat_ER_list(end)-mean(GCAS_GFAS_FluxSat_ER_list(1:end-1));
GCAS_GFAS_Mean_GPP_ER_2023anomaly=GCAS_GFAS_Mean_GPP_ER_list(end)-mean(GCAS_GFAS_Mean_GPP_ER_list(1:end-1));
GCAS_GFED_GOSIF_ER_2023anomaly=GCAS_GFED_GOSIF_ER_list(end)-mean(GCAS_GFED_GOSIF_ER_list(1:end-1));
GCAS_GFED_FluxSat_ER_2023anomaly=GCAS_GFED_FluxSat_ER_list(end)-mean(GCAS_GFED_FluxSat_ER_list(1:end-1));
GCAS_GFED_Mean_GPP_ER_2023anomaly=GCAS_GFED_Mean_GPP_ER_list(end)-mean(GCAS_GFED_Mean_GPP_ER_list(1:end-1));
GCAS_Mean_fire_GOSIF_ER_2023anomaly=GCAS_Mean_fire_GOSIF_ER_list(end)-mean(GCAS_Mean_fire_GOSIF_ER_list(1:end-1));
GCAS_Mean_fire_FluxSat_ER_2023anomaly=GCAS_Mean_fire_FluxSat_ER_list(end)-mean(GCAS_Mean_fire_FluxSat_ER_list(1:end-1));

GCAS_landflux_2023anomaly_std=std([GCAS_landflux_list(1,end)-mean(GCAS_landflux_list(1,1:end-1)),GCAS_landflux_list(2,end)-mean(GCAS_landflux_list(2,1:end-1)),...
    GCAS_landflux_list(3,end)-mean(GCAS_landflux_list(3,1:end-1)),GCAS_landflux_list(4,end)-mean(GCAS_landflux_list(4,1:end-1))]);
GCAS_GFAS_NEE_2023anomaly_std=sqrt(power((GFAS_fire_list(end)-mean(GFAS_fire_list(1:end-1)))*0.2,2)+power(GCAS_landflux_2023anomaly_std,2));
GCAS_GFED_NEE_2023anomaly_std=sqrt(power((GFED_fire_list(end)-mean(GFED_fire_list(1:end-1)))*0.2,2)+power(GCAS_landflux_2023anomaly_std,2));
GCAS_Mean_fire_NEE_2023anomaly_std=sqrt(power(Mean_fire_2023anomaly*0.2,2)+power(GCAS_landflux_2023anomaly_std,2));
GPP_2023anomaly_std=std(bootstrp(1000,@mean,[GOSIF_GPP_list(end)-mean(GOSIF_GPP_list(1:end-1)),FluxSat_GPP_list(end)-mean(FluxSat_GPP_list(1:end-1))]));
Fire_2023anomaly_std=Mean_fire_list(end)*0.2-mean(Mean_fire_list(1:end-1))*0.2;
GCAS_ER_2023anomaly_std=sqrt(power(Fire_2023anomaly_std,2)+power(GCAS_landflux_2023anomaly_std,2)+power(GPP_2023anomaly_std,2));

% GCAS 2023 ER martix
GCAS_ER_2023_martix=[GCAS_GFED_GOSIF_ER_2023,GCAS_GFED_Mean_GPP_ER_2023,GCAS_GFED_FluxSat_ER_2023;...
    GCAS_Mean_fire_GOSIF_ER_2023,GCAS_ER_2023,GCAS_Mean_fire_FluxSat_ER_2023;...
    GCAS_GFAS_GOSIF_ER_2023,GCAS_GFAS_Mean_GPP_ER_2023,GCAS_GFAS_FluxSat_ER_2023...
    ];
% GCAS 2023 ER anomaly martix
GCAS_ER_2023anomaly_martix=[GCAS_GFED_GOSIF_ER_2023anomaly,GCAS_GFED_Mean_GPP_ER_2023anomaly,GCAS_GFED_FluxSat_ER_2023anomaly;...
    GCAS_Mean_fire_GOSIF_ER_2023anomaly,GCAS_ER_2023anomaly,GCAS_Mean_fire_FluxSat_ER_2023anomaly;...
    GCAS_GFAS_GOSIF_ER_2023anomaly,GCAS_GFAS_Mean_GPP_ER_2023anomaly,GCAS_GFAS_FluxSat_ER_2023anomaly...
    ];
%% calculate landflux value (GCB2024-satellite)

for year=2015:2023

    % original fire
    GFAS_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\total_carbon\year\Fire_' num2str(year) '.tif']);
    GFED_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\total_carbon\year\Fire_' num2str(year) '.tif']);
    Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\year\Fire_' num2str(year) '.tif']);
    % case NEE
    CAMSS_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']);
    CMSF_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']);
    GCASv2_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']);
    GONGGA_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']);
    COLA_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\year\COLA_NEE_' num2str(year) '.tif']);
    NTFVAR_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\year\NTFVAR_NEE_' num2str(year) '.tif']);

    % original landflux
    CAMSS_landflux=CAMSS_mean_NEE+Mean_fire;
    CMSF_landflux=CMSF_mean_NEE+Mean_fire;
    GCASv2_landflux=GCASv2_mean_NEE+Mean_fire;
    GONGGA_landflux=GONGGA_mean_NEE+Mean_fire;
    COLA_landflux=COLA_mean_NEE+Mean_fire;
    NTFVAR_landflux=NTFVAR_mean_NEE+Mean_fire;

   % other case NEE
    CAMSS_GFED_NEE=CAMSS_landflux-GFED_fire;
    CAMSS_GFAS_NEE=CAMSS_landflux-GFAS_fire;
    CMSF_GFED_NEE=CMSF_landflux-GFED_fire;
    CMSF_GFAS_NEE=CMSF_landflux-GFAS_fire;
    GCASv2_GFED_NEE=GCASv2_landflux-GFED_fire;
    GCASv2_GFAS_NEE=GCASv2_landflux-GFAS_fire;
    GONGGA_GFED_NEE=GONGGA_landflux-GFED_fire;
    GONGGA_GFAS_NEE=GONGGA_landflux-GFAS_fire;
    COLA_GFED_NEE=COLA_landflux-GFED_fire;
    COLA_GFAS_NEE=COLA_landflux-GFAS_fire;
    NTFVAR_GFED_NEE=NTFVAR_landflux-GFED_fire;
    NTFVAR_GFAS_NEE=NTFVAR_landflux-GFAS_fire;

    % landflux list
    GCB_satellite_landflux_list(1,year-2014)=nansum(nansum(CAMSS_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(2,year-2014)=nansum(nansum(CMSF_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(3,year-2014)=nansum(nansum(GCASv2_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(4,year-2014)=nansum(nansum(GONGGA_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(5,year-2014)=nansum(nansum(COLA_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(6,year-2014)=nansum(nansum(NTFVAR_landflux.*area_grid/(10^15)));

    % GFAS-NEE list
    GCB_satellite_GFAS_NEE_list(1,year-2014)=nansum(nansum(CAMSS_GFAS_NEE.*area_grid/(10^15)));
    GCB_satellite_GFAS_NEE_list(2,year-2014)=nansum(nansum(CMSF_GFAS_NEE.*area_grid/(10^15)));
    GCB_satellite_GFAS_NEE_list(3,year-2014)=nansum(nansum(GCASv2_GFAS_NEE.*area_grid/(10^15)));
    GCB_satellite_GFAS_NEE_list(4,year-2014)=nansum(nansum(GONGGA_GFAS_NEE.*area_grid/(10^15)));
    GCB_satellite_GFAS_NEE_list(5,year-2014)=nansum(nansum(COLA_GFAS_NEE.*area_grid/(10^15)));
    GCB_satellite_GFAS_NEE_list(6,year-2014)=nansum(nansum(NTFVAR_GFAS_NEE.*area_grid/(10^15)));

    % GFED-NEE list
    GCB_satellite_GFED_NEE_list(1,year-2014)=nansum(nansum(CAMSS_GFED_NEE.*area_grid/(10^15)));
    GCB_satellite_GFED_NEE_list(2,year-2014)=nansum(nansum(CMSF_GFED_NEE.*area_grid/(10^15)));
    GCB_satellite_GFED_NEE_list(3,year-2014)=nansum(nansum(GCASv2_GFED_NEE.*area_grid/(10^15)));
    GCB_satellite_GFED_NEE_list(4,year-2014)=nansum(nansum(GONGGA_GFED_NEE.*area_grid/(10^15)));
    GCB_satellite_GFED_NEE_list(5,year-2014)=nansum(nansum(COLA_GFED_NEE.*area_grid/(10^15)));
    GCB_satellite_GFED_NEE_list(6,year-2014)=nansum(nansum(NTFVAR_GFED_NEE.*area_grid/(10^15)));
    
    % Meanfire-NEE list
    GCB_satellite_Mean_fire_NEE_list(1,year-2014)=nansum(nansum(CAMSS_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(2,year-2014)=nansum(nansum(CMSF_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(3,year-2014)=nansum(nansum(GCASv2_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(4,year-2014)=nansum(nansum(GONGGA_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(5,year-2014)=nansum(nansum(COLA_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(6,year-2014)=nansum(nansum(NTFVAR_mean_NEE.*area_grid/(10^15)));

end

% ER list
GCB_satellite_GFAS_GOSIF_ER_list=mean(GCB_satellite_GFAS_NEE_list)+GOSIF_GPP_list;
GCB_satellite_GFAS_FluxSat_ER_list=mean(GCB_satellite_GFAS_NEE_list)+FluxSat_GPP_list;
GCB_satellite_GFAS_Mean_GPP_ER_list=mean(GCB_satellite_GFAS_NEE_list)+Mean_GPP_list;
GCB_satellite_GFED_GOSIF_ER_list=mean(GCB_satellite_GFED_NEE_list)+GOSIF_GPP_list;
GCB_satellite_GFED_FluxSat_ER_list=mean(GCB_satellite_GFED_NEE_list)+FluxSat_GPP_list;
GCB_satellite_GFED_Mean_GPP_ER_list=mean(GCB_satellite_GFED_NEE_list)+Mean_GPP_list;
GCB_satellite_Mean_fire_GOSIF_ER_list=mean(GCB_satellite_Mean_fire_NEE_list)+GOSIF_GPP_list;
GCB_satellite_Mean_fire_FluxSat_ER_list=mean(GCB_satellite_Mean_fire_NEE_list)+FluxSat_GPP_list;
GCB_satellite_Mean_fire_Mean_GPP_ER_list=mean(GCB_satellite_Mean_fire_NEE_list)+Mean_GPP_list;

% 2023 Flux value
GCB_satellite_landflux_2023=mean(GCB_satellite_landflux_list(:,end));
GCB_satellite_Byrne_NEE_2023=GCB_satellite_landflux_2023-0.647;
GCB_satellite_GFAS_NEE_2023mean=mean(GCB_satellite_GFAS_NEE_list(:,end));
GCB_satellite_GFED_NEE_2023mean=mean(GCB_satellite_GFED_NEE_list(:,end));
GCB_satellite_Mean_fire_NEE_2023mean=mean(GCB_satellite_Mean_fire_NEE_list(:,end));
GCB_satellite_ER_2023=GCB_satellite_Mean_fire_Mean_GPP_ER_list(end);

GCB_satellite_GFAS_GOSIF_ER_2023=GCB_satellite_GFAS_GOSIF_ER_list(end);
GCB_satellite_GFAS_FluxSat_ER_2023=GCB_satellite_GFAS_FluxSat_ER_list(end);
GCB_satellite_GFAS_Mean_GPP_ER_2023=GCB_satellite_GFAS_Mean_GPP_ER_list(end);
GCB_satellite_GFED_GOSIF_ER_2023=GCB_satellite_GFED_GOSIF_ER_list(end);
GCB_satellite_GFED_FluxSat_ER_2023=GCB_satellite_GFED_FluxSat_ER_list(end);
GCB_satellite_GFED_Mean_GPP_ER_2023=GCB_satellite_GFED_Mean_GPP_ER_list(end);
GCB_satellite_Mean_fire_GOSIF_ER_2023=GCB_satellite_Mean_fire_GOSIF_ER_list(end);
GCB_satellite_Mean_fire_FluxSat_ER_2023=GCB_satellite_Mean_fire_FluxSat_ER_list(end);
GCB_satellite_Bynre_GOSIF_ER_2023=GCB_satellite_Byrne_NEE_2023+GOSIF_GPP_list(end);
GCB_satellite_Bynre_FluxSat_ER_2023=GCB_satellite_Byrne_NEE_2023+FluxSat_GPP_list(end);
GCB_satellite_Bynre_Mean_GPP_ER_2023=GCB_satellite_Byrne_NEE_2023+Mean_GPP_2023;


GCB_satellite_landflux_2023std=std(GCB_satellite_landflux_list(:,end));
GCB_satellite_GFAS_NEE_2023std=sqrt(power(GFAS_fire_list(end)*0.2,2)+power(GCB_satellite_landflux_2023std,2));
GCB_satellite_GFED_NEE_2023std=sqrt(power(GFED_fire_list(end)*0.2,2)+power(GCB_satellite_landflux_2023std,2));
GCB_satellite_Bynre_NEE_2023std=sqrt(power(0.647*0.2,2)+power(GCB_satellite_landflux_2023std,2));
GCB_satellite_Mean_fire_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(GCB_satellite_landflux_2023std,2));
GCB_satellite_ER_2023std=sqrt(power(Fire_2023std,2)+power(GCB_satellite_landflux_2023std,2)+power(GPP_2023std,2));

% 2023 Anomaly relative to 2015-2022
GCB_satellite_landflux_2023anomaly_mean=mean(GCB_satellite_landflux_list(:,end))-mean(mean(GCB_satellite_landflux_list(:,1:end-1)));
GCB_satellite_GFAS_NEE_2023anomaly_mean=mean(GCB_satellite_GFAS_NEE_list(:,end))-mean(mean(GCB_satellite_GFAS_NEE_list(:,1:end-1)));
GCB_satellite_GFED_NEE_2023anomaly_mean=mean(GCB_satellite_GFED_NEE_list(:,end))-mean(mean(GCB_satellite_GFED_NEE_list(:,1:end-1)));
GCB_satellite_Mean_fire_NEE_2023anomaly_mean=mean(GCB_satellite_Mean_fire_NEE_list(:,end))-mean(mean(GCB_satellite_Mean_fire_NEE_list(:,1:end-1)));
Mean_GPP_2023anomaly=Mean_GPP_list(end)-mean(Mean_GPP_list(1:end-1));
GCB_satellite_ER_2023anomaly=GCB_satellite_Mean_fire_Mean_GPP_ER_list(end)-mean(GCB_satellite_Mean_fire_Mean_GPP_ER_list(1:end-1));
Mean_fire_2023anomaly=Mean_fire_list(end)-mean(Mean_fire_list(1:end-1));

GCB_satellite_GFAS_GOSIF_ER_2023anomaly=GCB_satellite_GFAS_GOSIF_ER_list(end)-mean(GCB_satellite_GFAS_GOSIF_ER_list(1:end-1));
GCB_satellite_GFAS_FluxSat_ER_2023anomaly=GCB_satellite_GFAS_FluxSat_ER_list(end)-mean(GCB_satellite_GFAS_FluxSat_ER_list(1:end-1));
GCB_satellite_GFAS_Mean_GPP_ER_2023anomaly=GCB_satellite_GFAS_Mean_GPP_ER_list(end)-mean(GCB_satellite_GFAS_Mean_GPP_ER_list(1:end-1));
GCB_satellite_GFED_GOSIF_ER_2023anomaly=GCB_satellite_GFED_GOSIF_ER_list(end)-mean(GCB_satellite_GFED_GOSIF_ER_list(1:end-1));
GCB_satellite_GFED_FluxSat_ER_2023anomaly=GCB_satellite_GFED_FluxSat_ER_list(end)-mean(GCB_satellite_GFED_FluxSat_ER_list(1:end-1));
GCB_satellite_GFED_Mean_GPP_ER_2023anomaly=GCB_satellite_GFED_Mean_GPP_ER_list(end)-mean(GCB_satellite_GFED_Mean_GPP_ER_list(1:end-1));
GCB_satellite_Mean_fire_GOSIF_ER_2023anomaly=GCB_satellite_Mean_fire_GOSIF_ER_list(end)-mean(GCB_satellite_Mean_fire_GOSIF_ER_list(1:end-1));
GCB_satellite_Mean_fire_FluxSat_ER_2023anomaly=GCB_satellite_Mean_fire_FluxSat_ER_list(end)-mean(GCB_satellite_Mean_fire_FluxSat_ER_list(1:end-1));

GCB_satellite_landflux_2023anomaly_std=std([GCB_satellite_landflux_list(1,end)-mean(GCB_satellite_landflux_list(1,1:end-1)),GCB_satellite_landflux_list(2,end)-mean(GCB_satellite_landflux_list(2,1:end-1)),...
    GCB_satellite_landflux_list(3,end)-mean(GCB_satellite_landflux_list(3,1:end-1)),GCB_satellite_landflux_list(4,end)-mean(GCB_satellite_landflux_list(4,1:end-1))]);
GCB_satellite_GFAS_NEE_2023anomaly_std=sqrt(power((GFAS_fire_list(end)-mean(GFAS_fire_list(1:end-1)))*0.2,2)+power(GCB_satellite_landflux_2023anomaly_std,2));
GCB_satellite_GFED_NEE_2023anomaly_std=sqrt(power((GFED_fire_list(end)-mean(GFED_fire_list(1:end-1)))*0.2,2)+power(GCB_satellite_landflux_2023anomaly_std,2));
GCB_satellite_Mean_fire_NEE_2023anomaly_std=sqrt(power(Mean_fire_2023anomaly*0.2,2)+power(GCB_satellite_landflux_2023anomaly_std,2));
GPP_2023anomaly_std=std(bootstrp(1000,@mean,[GOSIF_GPP_list(end)-mean(GOSIF_GPP_list(1:end-1)),FluxSat_GPP_list(end)-mean(FluxSat_GPP_list(1:end-1))]));
Fire_2023anomaly_std=Mean_fire_list(end)*0.2-mean(Mean_fire_list(1:end-1))*0.2;
GCB_satellite_ER_2023anomaly_std=sqrt(power(Fire_2023anomaly_std,2)+power(GCB_satellite_landflux_2023anomaly_std,2)+power(GPP_2023anomaly_std,2));


% GCB_satellite 2023 ER  martix
GCB_satellite_ER_2023_martix=[
    GCB_satellite_GFED_GOSIF_ER_2023,GCB_satellite_GFED_Mean_GPP_ER_2023,GCB_satellite_GFED_FluxSat_ER_2023;...
    GCB_satellite_Mean_fire_GOSIF_ER_2023,GCB_satellite_ER_2023,GCB_satellite_Mean_fire_FluxSat_ER_2023;...
    GCB_satellite_GFAS_GOSIF_ER_2023,GCB_satellite_GFAS_Mean_GPP_ER_2023,GCB_satellite_GFAS_FluxSat_ER_2023...
    ];
% GCB_satellite 2023 ER anomaly  martix
GCB_satellite_ER_2023anomaly_martix=[
    GCB_satellite_GFED_GOSIF_ER_2023anomaly,GCB_satellite_GFED_Mean_GPP_ER_2023anomaly,GCB_satellite_GFED_FluxSat_ER_2023anomaly;
    GCB_satellite_Mean_fire_GOSIF_ER_2023anomaly,GCB_satellite_ER_2023anomaly,GCB_satellite_Mean_fire_FluxSat_ER_2023anomaly;...
    GCB_satellite_GFAS_GOSIF_ER_2023anomaly,GCB_satellite_GFAS_Mean_GPP_ER_2023anomaly,GCB_satellite_GFAS_FluxSat_ER_2023anomaly...
    ];
%% total Fire (CO2+CO)

for year=2015:2023

    GFED_total_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\year\GFED_' num2str(year) '.tif']);
    GFAS_total_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\CO+CH4\year\Fire_CO+CH4_' num2str(year) '.tif'])+importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\year\GFAS_' num2str(year) '.tif']);
    Mean_total_fire=(GFED_total_fire+GFAS_total_fire)/2;

    GFED_total_fire_list(year-2014)=nansum(nansum(GFED_total_fire.*area_grid/(10^15)));
    GFAS_total_fire_list(year-2014)=nansum(nansum(GFAS_total_fire.*area_grid/(10^15)));
    Mean_total_fire_list(year-2014)=nansum(nansum(Mean_total_fire.*area_grid/(10^15)));
end

%%
total_landflux_2023std=std([GCB_satellite_landflux_list(:,end);GCAS_landflux_list(:,end)]);
total_NEE_2023=mean([GCB_satellite_Mean_fire_NEE_list(:,end);GCAS_Mean_fire_NEE_list(:,end)]);
total_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(total_landflux_2023std,2));


total_landflux_2023anomaly_std=std([GCAS_landflux_list(1,end)-mean(GCAS_landflux_list(1,1:end-1)),GCAS_landflux_list(2,end)-mean(GCAS_landflux_list(2,1:end-1)),...
    GCAS_landflux_list(3,end)-mean(GCAS_landflux_list(3,1:end-1)),GCAS_landflux_list(4,end)-mean(GCAS_landflux_list(4,1:end-1)),...
    GCB_satellite_landflux_list(1,end)-mean(GCB_satellite_landflux_list(1,1:end-1)),GCB_satellite_landflux_list(2,end)-mean(GCB_satellite_landflux_list(2,1:end-1)),...
    GCB_satellite_landflux_list(3,end)-mean(GCB_satellite_landflux_list(3,1:end-1)),GCB_satellite_landflux_list(4,end)-mean(GCB_satellite_landflux_list(4,1:end-1))]);
total_NEE_2023anomaly_mean=mean([GCAS_Mean_fire_NEE_list(:,end);GCB_satellite_Mean_fire_NEE_list(:,end)])-mean(mean([GCAS_Mean_fire_NEE_list(:,1:end-1);GCB_satellite_Mean_fire_NEE_list(:,1:end-1)]));

total_NEE_2023anomalystd=sqrt(power(Fire_2023anomaly_std(end)*0.2,2)+power(total_landflux_2023anomaly_std,2));

% contribution rate
-GCAS_Mean_fire_NEE_2023anomaly_mean/Mean_fire_2023(end)
% -GCAS_Mean_fire_NEE_2023anomaly_mean/0.647
-GCB_satellite_Mean_fire_NEE_2023anomaly_mean/Mean_fire_2023(end)
% -GCB_satellite_Mean_fire_NEE_2023anomaly_mean/0.647


Mean_GPP_2023anomaly/GCAS_Mean_fire_NEE_2023anomaly_mean
Mean_GPP_2023anomaly/GCB_satellite_Mean_fire_NEE_2023anomaly_mean

%%

ff=figure;
set(gcf,'unit','pixels','position',[1000,892,1035,446]);

% part1---------------------------------------------------------------
a=axes('Position',[0.0574676328502415 0.344755784061697 0.421355959016828 0.619315570762208]);
b1=bar(1-0.25,-GCAS_landflux_2023,'BarWidth',0.5,'FaceColor',[93,84,113]/255);hold on 
errorbar(1-0.25,-GCAS_landflux_2023,GCAS_landflux_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b2=bar(1+0.25,-GCB_satellite_landflux_2023,'BarWidth',0.5,'FaceColor','k');hold on
errorbar(1+0.25,-GCB_satellite_landflux_2023,GCB_satellite_landflux_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

b3=bar(3-0.25,-GCAS_GFED_NEE_2023mean,'BarWidth',0.5,'FaceColor',[197,223,248]/255);hold on
errorbar(3-0.25,-GCAS_GFED_NEE_2023mean,GCAS_GFED_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b4=bar(3+0.25,-GCAS_GFAS_NEE_2023mean,'BarWidth',0.5,'FaceColor',[161,191,224]/255);hold on
errorbar(3+0.25,-GCAS_GFAS_NEE_2023mean,GCAS_GFAS_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b5=bar(4-0.25,-GCAS_Mean_fire_NEE_2023mean,'BarWidth',0.5,'FaceColor',[120,149,203]/255);hold on
errorbar(4-0.25,-GCAS_Mean_fire_NEE_2023mean,GCAS_Mean_fire_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b6=bar(4+0.25,-GCB_satellite_GFED_NEE_2023mean,'BarWidth',0.5,'FaceColor',[254,227,229]/255);hold on
errorbar(4+0.25,-GCB_satellite_GFED_NEE_2023mean,GCB_satellite_GFED_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b7=bar(5-0.25,-GCB_satellite_GFAS_NEE_2023mean,'BarWidth',0.5,'FaceColor',[236,169,207]/255);hold on
errorbar(5-0.25,-GCB_satellite_GFAS_NEE_2023mean,GCB_satellite_GFAS_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b8=bar(5+0.25,-GCB_satellite_Mean_fire_NEE_2023mean,'BarWidth',0.5,'FaceColor',[184,128,197]/255);hold on
errorbar(5+0.25,-GCB_satellite_Mean_fire_NEE_2023mean,GCB_satellite_Mean_fire_NEE_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

b9=bar(7-0.25,Mean_GPP_2023,'BarWidth',0.5,'FaceColor',[17,119,51]/255);hold on
errorbar(7-0.25,Mean_GPP_2023,GPP_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b10=bar(8.5-0.25,GCAS_ER_2023,'BarWidth',0.5,'FaceColor',[141,124,106]/255);hold on
errorbar(8.5-0.25,GCAS_ER_2023,GCAS_ER_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b11=bar(8.5+0.25,GCB_satellite_ER_2023,'BarWidth',0.5,'FaceColor',[200,184,147]/255);hold on
errorbar(8.5+0.25,GCAS_ER_2023,GCB_satellite_ER_2023std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')


b12=bar(10.25-0.25,GFED_total_fire_list(end),'BarWidth',0.5,'FaceColor',[252,180,168]/255);hold on
errorbar(10.25-0.25,GFED_total_fire_list(end),GFED_total_fire_list(end)*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b13=bar(10.25+0.25,GFAS_total_fire_list(end),'BarWidth',0.5,'FaceColor',[247,122,130]/255);hold on
errorbar(10.25+0.25,GFAS_total_fire_list(end),GFAS_total_fire_list(end)*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b14=bar(11.25-0.25,Mean_total_fire_list(end),'BarWidth',0.5,'FaceColor',[219,93,107]/255);hold on
errorbar(11.25-0.25,Mean_total_fire_list(end),Mean_total_fire_list(end)*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
b15=bar(11.25+0.25,0.647,'BarWidth',0.5,'FaceColor',[200,55,86]/255);hold on
errorbar(11.25+0.25,0.647,0.727-0.647,0.647-0.57, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

xlim([0,12.25])
set(gca,'XTick',[1,4,6.75,8.5,10.75],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'NBP','NEP','GPP','TER','Fire'},'FontName','Arial','fontsize',14)
set(gca,'YTick',[-1:1:4],'FontName','Arial','fontsize',14)
% ylabel('PgC','FontName','Arial','FontSize',14);
truncAxis('Y',[1.3,2.6])
ylabel('PgC','FontName','Arial','FontSize',14,'position',[-0.9164 2.1370]);

lgd=legend([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15],{'NBP (GCAS-extra)','NBP (GCB2024-satellite)','NEP (GCAS-extra NBP + GFED fire)',...
    'NEP (GCAS-extra NBP + GFAS fire)','NEP (GCAS-extra NBP + MGG fire)','NEP (GCB2024-satellite NBP + GFED fire)',...
    'NEP (GCB2024-satellite NBP + GFAS fire)','NEP (GCB2024-satellite NBP + MGG fire)',...
    'GPP','TER (GCAS-extra)','TER (GCB2024-satellite)','Fire (GFED)','Fire (GFAS)','Mean of GFED and GFAS (MGG)','Fire (Byrne)'}...
    ,'NumColumns',3,'FontName','Arial','FontSize',13,'Box','off','Location','southoutside')
% lgd.ItemTokenSize = [9 9];
lgd.Position=[0.00121963808302931 0.020206803641003 0.997101422208518 0.209511562300832]

text('string','a','Units','normalized','position',[-0.126824003590029 1.00478120967781 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part2---------------------------------------------------------------
b=axes('Position',[0.563185103785104 0.344755784061697 0.421355959016828 0.619315570762208]);
bar(1-0.25,-GCAS_landflux_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[93,84,113]/255);hold on 
errorbar(1-0.25,-GCAS_landflux_2023anomaly_mean,GCAS_landflux_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(1+0.25,-GCB_satellite_landflux_2023anomaly_mean,'BarWidth',0.5,'FaceColor','k');hold on
errorbar(1+0.25,-GCB_satellite_landflux_2023anomaly_mean,GCB_satellite_landflux_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

bar(3-0.25,-GCAS_GFAS_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[197,223,248]/255);hold on
errorbar(3-0.25,-GCAS_GFAS_NEE_2023anomaly_mean,GCAS_GFAS_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(3+0.25,-GCAS_GFED_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[161,191,224]/255);hold on
errorbar(3+0.25,-GCAS_GFED_NEE_2023anomaly_mean,GCAS_GFED_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(4-0.25,-GCAS_Mean_fire_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[120,149,203]/255);hold on
errorbar(4-0.25,-GCAS_Mean_fire_NEE_2023anomaly_mean,GCAS_Mean_fire_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

bar(4+0.25,-GCB_satellite_GFAS_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[254,227,229]/255);hold on
errorbar(4+0.25,-GCB_satellite_GFAS_NEE_2023anomaly_mean,GCB_satellite_GFAS_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(5-0.25,-GCB_satellite_GFED_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[236,169,207]/255);hold on
errorbar(5-0.25,-GCB_satellite_GFED_NEE_2023anomaly_mean,GCB_satellite_GFED_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(5+0.25,-GCB_satellite_Mean_fire_NEE_2023anomaly_mean,'BarWidth',0.5,'FaceColor',[184,128,197]/255);hold on
errorbar(5+0.25,-GCB_satellite_Mean_fire_NEE_2023anomaly_mean,GCB_satellite_Mean_fire_NEE_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

bar(7-0.25,Mean_GPP_2023anomaly,'BarWidth',0.5,'FaceColor',[17,119,51]/255);hold on
errorbar(7-0.25,Mean_GPP_2023anomaly,GPP_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(8.5-0.25,GCAS_ER_2023anomaly,'BarWidth',0.5,'FaceColor',[141,124,106]/255);hold on
errorbar(8.5-0.25,GCAS_ER_2023anomaly,GCAS_ER_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(8.5+0.25,GCB_satellite_ER_2023anomaly,'BarWidth',0.5,'FaceColor',[200,184,147]/255);hold on
errorbar(8.5+0.25,GCB_satellite_ER_2023anomaly,GCB_satellite_ER_2023anomaly_std, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')


bar(10.5-0.25,GFED_total_fire_list(end)-mean(GFED_total_fire_list(1:end-1)),'BarWidth',0.5,'FaceColor',[252,180,168]/255);hold on
errorbar(10.5-0.25,GFED_total_fire_list(end)-mean(GFED_total_fire_list(1:end-1)),(GFED_total_fire_list(end)-mean(GFED_total_fire_list(1:end-1)))*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(10.5+0.25,GFAS_total_fire_list(end)-mean(GFAS_total_fire_list(1:end-1)),'BarWidth',0.5,'FaceColor',[247,122,130]/255);hold on
errorbar(10.5+0.25,GFAS_total_fire_list(end)-mean(GFAS_total_fire_list(1:end-1)),(GFAS_total_fire_list(end)-mean(GFAS_total_fire_list(1:end-1)))*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')
bar(11.5-0.25,Mean_total_fire_list(end)-mean(Mean_total_fire_list(1:end-1)),'BarWidth',0.5,'FaceColor',[219,93,107]/255);hold on
errorbar(11.5-0.25,Mean_total_fire_list(end)-mean(Mean_total_fire_list(1:end-1)),(Mean_total_fire_list(end)-mean(Mean_total_fire_list(1:end-1)))*0.2, ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

xlim([0,12])
set(gca,'XTick',[1,4,6.75,8.5,10.75],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'NBP','NEP','GPP','TER','Fire'},'FontName','Arial','fontsize',14)
ylabel('\DeltaPgC','FontName','Arial','FontSize',14);
ylim([-0.7,1])
text('string','b','Units','normalized','position',[-0.126824003590029 1.00478120967781 0],'FontName','Arial','FontSize',20,'fontweight','bold')

result=['E:\phd_file\Boreal_North_America\Result\V9\Fire_NEE_bar.png']
% print(result,ff,'-r600','-dpng' );
%%
ff=figure;
set(gcf,'unit','pixels','position',[1000,725,726,513]);
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% part1---------------------------------------------------------------
nexttile
% c=axes('Position',[0.0610618357487918 0.03257648387945 0.166638647342996,0.338350724637681]);
h1 = heatmap2(GCAS_ER_2023_martix,{'GOSIF','Mean','FluxSat'},{'GFED','MGG','GFAS'},'%0.2f','GridLines', '-');
caxis([2.6,3.6]);
ax=gca;
ax.FontSize=12;
ax.FontName='Arial';
colormap(nclCM(463))
c1=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[2.6:0.5:3.6]);
set(get(c1,'ylabel'),'string','TER (PgC)','fontsize',14);
text('string','a','Units','normalized','position',[-0.244471062413558 1.02495840161774 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% d=axes('Position',[0.313264734299516 0.03257648387945 0.166638647342996,0.338350724637681]);
nexttile
h2 = heatmap2(GCB_satellite_ER_2023_martix,{'GOSIF','Mean','FluxSat'},{'GFED','MGG','GFAS'},'%0.2f','GridLines', '-');
caxis([2.6,3.6]);
ax=gca;
ax.FontSize=12;
ax.FontName='Arial';
colormap(nclCM(463))
c2=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[2.6:0.5:3.6]);
set(get(c2,'ylabel'),'string','TER (PgC)','fontsize',14);
text('string','b','Units','normalized','position',[-0.244471062413558 1.02495840161774 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% e=axes('Position',[0.564907246376811 0.03257648387945 0.166638647342996,0.338350724637681]);
e=nexttile
h3 = heatmap2(GCAS_ER_2023anomaly_martix,{'GOSIF','Mean','FluxSat'},{'GFED','MGG','GFAS'},'%0.2f','GridLines', '-');
caxis([-0.5,-0.1]);
ax=gca;
ax.FontSize=12;
ax.FontName='Arial';
colormap(e,flipud(nclCM(463)))
c3=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.5:0.2:-0.1]);
set(get(c3,'ylabel'),'string','TER anomalies (PgC)','fontsize',14);
text('string','c','Units','normalized','position',[-0.244471062413558 1.02495840161774 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% f=axes('Position',[0.807467632850241 0.03257648387945 0.166638647342996,0.338350724637681]);
f=nexttile
h4 = heatmap2(GCB_satellite_ER_2023anomaly_martix,{'GOSIF','Mean','FluxSat'},{'GFED','MGG','GFAS'},'%0.2f','GridLines', '-');
caxis([-0.5,-0.1]);
ax=gca;
ax.FontSize=12;
ax.FontName='Arial';
colormap(f,flipud(nclCM(463)))
c4=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.5:0.2:-0.1]);
set(get(c4,'ylabel'),'string','TER anomalies (PgC)','fontsize',14);
text('string','d','Units','normalized','position',[-0.244471062413558 1.02495840161774 0],'FontName','Arial','FontSize',20,'fontweight','bold')
colormap(e,flipud(nclCM(463)))

% set(gcf,'unit','centimeters','position',[14.605000000000002,6.2865,34.798,20.277666666666672]);
result=['E:\phd_file\Boreal_North_America\Result\V9\Fire_NEE_bar2.png']
% print(result,ff,'-r600','-dpng' );

