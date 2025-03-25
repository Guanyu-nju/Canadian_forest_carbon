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
GCAS_Mean_fire_Mean_GPP_ER_list=mean(GCAS_Mean_fire_NEE_list)+Mean_GPP_list;

% 2023 Flux value
GCAS_landflux_2023=mean(GCAS_landflux_list(:,end));
GCAS_Mean_fire_NEE_2023mean=mean(GCAS_Mean_fire_NEE_list(:,end));
GCAS_ER_2023=GCAS_Mean_fire_Mean_GPP_ER_list(end);

GCAS_landflux_2023std=std(GCAS_landflux_list(:,end));
GCAS_Mean_fire_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(GCAS_landflux_2023std,2));
GPP_2023std=std(bootstrp(1000,@mean,[GOSIF_GPP_list(end),FluxSat_GPP_list(end)]));
Fire_2023std=Mean_fire_list(end)*0.2;
GCAS_ER_2023std=sqrt(power(Fire_2023std,2)+power(GCAS_landflux_2023std,2)+power(GPP_2023std,2));

% 2023 Anomaly relative to 2015-2022
GCAS_landflux_2023anomaly_mean=mean(GCAS_landflux_list(:,end))-mean(mean(GCAS_landflux_list(:,1:end-1)));
GCAS_Mean_fire_NEE_2023anomaly_mean=mean(GCAS_Mean_fire_NEE_list(:,end))-mean(mean(GCAS_Mean_fire_NEE_list(:,1:end-1)));
Mean_GPP_2023anomaly=Mean_GPP_list(end)-mean(Mean_GPP_list(1:end-1));
GCAS_ER_2023anomaly=GCAS_Mean_fire_Mean_GPP_ER_list(end)-mean(GCAS_Mean_fire_Mean_GPP_ER_list(1:end-1));
Mean_fire_2023anomaly=Mean_fire_list(end)-mean(Mean_fire_list(1:end-1));


GCAS_landflux_2023anomaly_std=std([GCAS_landflux_list(1,end)-mean(GCAS_landflux_list(1,1:end-1)),GCAS_landflux_list(2,end)-mean(GCAS_landflux_list(2,1:end-1)),...
    GCAS_landflux_list(3,end)-mean(GCAS_landflux_list(3,1:end-1)),GCAS_landflux_list(4,end)-mean(GCAS_landflux_list(4,1:end-1))]);
GCAS_Mean_fire_NEE_2023anomaly_std=sqrt(power(Mean_fire_2023anomaly*0.2,2)+power(GCAS_landflux_2023anomaly_std,2));
GPP_2023anomaly_std=std(bootstrp(1000,@mean,[GOSIF_GPP_list(end)-mean(GOSIF_GPP_list(1:end-1)),FluxSat_GPP_list(end)-mean(FluxSat_GPP_list(1:end-1))]));
Fire_2023anomaly_std=Mean_fire_list(end)*0.2-mean(Mean_fire_list(1:end-1))*0.2;
GCAS_ER_2023anomaly_std=sqrt(power(Fire_2023anomaly_std,2)+power(GCAS_landflux_2023anomaly_std,2)+power(GPP_2023anomaly_std,2));

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

    % landflux list
    GCB_satellite_landflux_list(1,year-2014)=nansum(nansum(CAMSS_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(2,year-2014)=nansum(nansum(CMSF_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(3,year-2014)=nansum(nansum(GCASv2_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(4,year-2014)=nansum(nansum(GONGGA_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(5,year-2014)=nansum(nansum(COLA_landflux.*area_grid/(10^15)));
    GCB_satellite_landflux_list(6,year-2014)=nansum(nansum(NTFVAR_landflux.*area_grid/(10^15)));

    % Meanfire-NEE list
    GCB_satellite_Mean_fire_NEE_list(1,year-2014)=nansum(nansum(CAMSS_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(2,year-2014)=nansum(nansum(CMSF_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(3,year-2014)=nansum(nansum(GCASv2_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(4,year-2014)=nansum(nansum(GONGGA_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(5,year-2014)=nansum(nansum(COLA_mean_NEE.*area_grid/(10^15)));
    GCB_satellite_Mean_fire_NEE_list(6,year-2014)=nansum(nansum(NTFVAR_mean_NEE.*area_grid/(10^15)));

end

% ER list
GCB_satellite_Mean_fire_Mean_GPP_ER_list=mean(GCB_satellite_Mean_fire_NEE_list)+Mean_GPP_list;

% 2023 Flux value
GCB_satellite_landflux_2023=mean(GCB_satellite_landflux_list(:,end));
GCB_satellite_Mean_fire_NEE_2023mean=mean(GCB_satellite_Mean_fire_NEE_list(:,end));
GCB_satellite_ER_2023=GCB_satellite_Mean_fire_Mean_GPP_ER_list(end);

GCB_satellite_landflux_2023std=std(GCB_satellite_landflux_list(:,end));
GCB_satellite_Mean_fire_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(GCB_satellite_landflux_2023std,2));
GCB_satellite_ER_2023std=sqrt(power(Fire_2023std,2)+power(GCB_satellite_landflux_2023std,2)+power(GPP_2023std,2));

% 2023 Anomaly relative to 2015-2022
GCB_satellite_landflux_2023anomaly_mean=mean(GCB_satellite_landflux_list(:,end))-mean(mean(GCB_satellite_landflux_list(:,1:end-1)));
GCB_satellite_Mean_fire_NEE_2023anomaly_mean=mean(GCB_satellite_Mean_fire_NEE_list(:,end))-mean(mean(GCB_satellite_Mean_fire_NEE_list(:,1:end-1)));
GCB_satellite_ER_2023anomaly=GCB_satellite_Mean_fire_Mean_GPP_ER_list(end)-mean(GCB_satellite_Mean_fire_Mean_GPP_ER_list(1:end-1));

GCB_satellite_landflux_2023anomaly_std=std([GCB_satellite_landflux_list(1,end)-mean(GCB_satellite_landflux_list(1,1:end-1)),GCB_satellite_landflux_list(2,end)-mean(GCB_satellite_landflux_list(2,1:end-1)),...
    GCB_satellite_landflux_list(3,end)-mean(GCB_satellite_landflux_list(3,1:end-1)),GCB_satellite_landflux_list(4,end)-mean(GCB_satellite_landflux_list(4,1:end-1))]);
GCB_satellite_Mean_fire_NEE_2023anomaly_std=sqrt(power(Mean_fire_2023anomaly*0.2,2)+power(GCB_satellite_landflux_2023anomaly_std,2));
GCB_satellite_ER_2023anomaly_std=sqrt(power(Fire_2023anomaly_std,2)+power(GCB_satellite_landflux_2023anomaly_std,2)+power(GPP_2023anomaly_std,2));
%% calculate landflux value (GCB2024-obspack)

for year=2015:2023

    % original fire
    Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\year\Fire_' num2str(year) '.tif']);

    % case NEE
    A_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\year\CAMS_NEE_' num2str(year) '.tif']);
    B_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\year\CarboScope_NEE_' num2str(year) '.tif']);
    C_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\year\CTE_NEE_' num2str(year) '.tif']);
    D_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\year\CT_NEE_' num2str(year) '.tif']);
    E_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\year\IAPCAS_NEE_' num2str(year) '.tif']);
    F_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\year\MIROC-ACTM_NEE_' num2str(year) '.tif']);
    G_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\year\NISMON_NEE_' num2str(year) '.tif']);
    H_mean_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\year\UoE_NEE_' num2str(year) '.tif']);

    % original landflux
    A_landflux=A_mean_NEE+Mean_fire;
    B_landflux=B_mean_NEE+Mean_fire;
    C_landflux=C_mean_NEE+Mean_fire;
    D_landflux=D_mean_NEE+Mean_fire;
    E_landflux=E_mean_NEE+Mean_fire;
    F_landflux=F_mean_NEE+Mean_fire;
    G_landflux=G_mean_NEE+Mean_fire;
    H_landflux=H_mean_NEE+Mean_fire;


    % landflux list
    GCB_obspack_landflux_list(1,year-2014)=nansum(nansum(A_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(2,year-2014)=nansum(nansum(B_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(3,year-2014)=nansum(nansum(C_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(4,year-2014)=nansum(nansum(D_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(5,year-2014)=nansum(nansum(E_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(6,year-2014)=nansum(nansum(F_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(7,year-2014)=nansum(nansum(G_landflux.*area_grid/(10^15)));
    GCB_obspack_landflux_list(8,year-2014)=nansum(nansum(H_landflux.*area_grid/(10^15)));


    % Meanfire-NEE list
    GCB_obspack_Mean_fire_NEE_list(1,year-2014)=nansum(nansum(A_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(2,year-2014)=nansum(nansum(B_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(3,year-2014)=nansum(nansum(C_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(4,year-2014)=nansum(nansum(D_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(5,year-2014)=nansum(nansum(E_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(6,year-2014)=nansum(nansum(F_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(7,year-2014)=nansum(nansum(G_mean_NEE.*area_grid/(10^15)));
    GCB_obspack_Mean_fire_NEE_list(8,year-2014)=nansum(nansum(H_mean_NEE.*area_grid/(10^15)));


end

% ER list
GCB_obspack_Mean_fire_Mean_GPP_ER_list=mean(GCB_obspack_Mean_fire_NEE_list)+Mean_GPP_list;

% 2023 Flux value
GCB_obspack_landflux_2023=mean(GCB_obspack_landflux_list(:,end));
GCB_obspack_Mean_fire_NEE_2023mean=mean(GCB_obspack_Mean_fire_NEE_list(:,end));
GCB_obspack_ER_2023=GCB_obspack_Mean_fire_Mean_GPP_ER_list(end);

GCB_obspack_landflux_2023std=std(GCB_obspack_landflux_list(:,end));
GCB_obspack_Mean_fire_NEE_2023std=sqrt(power(Mean_fire_list(end)*0.2,2)+power(GCB_obspack_landflux_2023std,2));
GCB_obspack_ER_2023std=sqrt(power(Fire_2023std,2)+power(GCB_obspack_landflux_2023std,2)+power(GPP_2023std,2));

% 2023 Anomaly relative to 2015-2022
GCB_obspack_landflux_2023anomaly_mean=mean(GCB_obspack_landflux_list(:,end))-mean(mean(GCB_obspack_landflux_list(:,1:end-1)));
GCB_obspack_Mean_fire_NEE_2023anomaly_mean=mean(GCB_obspack_Mean_fire_NEE_list(:,end))-mean(mean(GCB_obspack_Mean_fire_NEE_list(:,1:end-1)));
GCB_obspack_ER_2023anomaly=GCB_obspack_Mean_fire_Mean_GPP_ER_list(end)-mean(GCB_obspack_Mean_fire_Mean_GPP_ER_list(1:end-1));

GCB_obspack_landflux_2023anomaly_std=std([GCB_obspack_landflux_list(1,end)-mean(GCB_obspack_landflux_list(1,1:end-1)),GCB_obspack_landflux_list(2,end)-mean(GCB_obspack_landflux_list(2,1:end-1)),...
    GCB_obspack_landflux_list(3,end)-mean(GCB_obspack_landflux_list(3,1:end-1)),GCB_obspack_landflux_list(4,end)-mean(GCB_obspack_landflux_list(4,1:end-1))]);
GCB_obspack_Mean_fire_NEE_2023anomaly_std=sqrt(power(Mean_fire_2023anomaly*0.2,2)+power(GCB_obspack_landflux_2023anomaly_std,2));
GCB_obspack_ER_2023anomaly_std=sqrt(power(Fire_2023anomaly_std,2)+power(GCB_obspack_landflux_2023anomaly_std,2)+power(GPP_2023anomaly_std,2));
%% calculate mean value (GCB2024-obspack)

for year=2015:2022

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\year\NEE_' num2str(year) '.tif']);
    GCB_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\year\ER_' num2str(year) '.tif']);
    GCB_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);
GCB_ER_mean=nanmean(GCB_ER_mean,3);



GCB_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\year\NEE_2023.tif").*pixel_mask;
GCB_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\year\ER_2023.tif").*pixel_mask;


GCB_NEP_2023Anomaly=-GCB_NEE2023+GCB_NEE_mean;
GCB_ER_2023Anomaly=GCB_ER2023-GCB_ER_mean;
%% GCB2024

% 时间序列均值
for month=1:12


    for year=2015:2022

        Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Fire_list(year-2014,month)=nansum(nansum(Mean_fire.*pixel_mask.*area_grid/(10^15)));
        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GCB_NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));


        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year) '_' num2str(month) '.tif']);
        E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year) '_' num2str(month) '.tif']);
        G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year) '_' num2str(month) '.tif']);
        H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year) '_' num2str(month) '.tif']);

        % original landflux
        A_landflux=A_NEE+Mean_fire;
        B_landflux=B_NEE+Mean_fire;
        C_landflux=C_NEE+Mean_fire;
        D_landflux=D_NEE+Mean_fire;
        E_landflux=E_NEE+Mean_fire;
        F_landflux=F_NEE+Mean_fire;
        G_landflux=G_NEE+Mean_fire;
        H_landflux=H_NEE+Mean_fire;


        GCB_landflux_list(year-2014,month,1)=nansum(nansum(A_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,2)=nansum(nansum(B_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,3)=nansum(nansum(C_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,4)=nansum(nansum(D_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,5)=nansum(nansum(E_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,6)=nansum(nansum(F_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,7)=nansum(nansum(G_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,8)=nansum(nansum(H_landflux.*area_grid/(10^15)));


        GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
        GPP_list(year-2014,month)=nansum(nansum(GPP_mean_temp.*pixel_mask.*area_grid/(10^15)));

        GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
        Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

        GPP_std_list(year-2014,month,1)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
        GPP_std_list(year-2014,month,2)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));

        ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']);
        GCB_ER_list(year-2014,month)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    end

end

% 2023年数据
for month=1:12


    year=2023;
    Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    GCB_NEE_2023(month)=nansum(nansum(NNE_2023_temp.*pixel_mask.*area_grid/(10^15)));
    Fire_2023(month)=nansum(nansum(Mean_fire.*pixel_mask.*area_grid/(10^15)));



    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year) '_' num2str(month) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year) '_' num2str(month) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year) '_' num2str(month) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year) '_' num2str(month) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year) '_' num2str(month) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year) '_' num2str(month) '.tif']);
    G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year) '_' num2str(month) '.tif']);
    H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year) '_' num2str(month) '.tif']);

    % original landflux
    A_landflux=A_NEE+Mean_fire;
    B_landflux=B_NEE+Mean_fire;
    C_landflux=C_NEE+Mean_fire;
    D_landflux=D_NEE+Mean_fire;
    E_landflux=E_NEE+Mean_fire;
    F_landflux=F_NEE+Mean_fire;
    G_landflux=G_NEE+Mean_fire;
    H_landflux=H_NEE+Mean_fire;


    GCB_landflux_2023_list(1,month)=nansum(nansum(A_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(2,month)=nansum(nansum(B_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(3,month)=nansum(nansum(C_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(4,month)=nansum(nansum(D_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(5,month)=nansum(nansum(E_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(6,month)=nansum(nansum(F_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(7,month)=nansum(nansum(G_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(8,month)=nansum(nansum(H_landflux.*pixel_mask.*area_grid/(10^15)));


    GPP_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
    GPP_2023(month)=nansum(nansum(GPP_2023_temp.*pixel_mask.*area_grid/(10^15)));

    GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
    Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

    GPP_2023_list(1,month)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
    GPP_2023_list(2,month)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));


    ER_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']);
    GCB_ER_2023(month)=nansum(nansum(ER_2023_temp.*pixel_mask.*area_grid/(10^15)));


end

% 2023 flux anomaly
GCB_NEP_anomaly_mean=-GCB_NEE_2023+nanmean(GCB_NEE_list);

for i=1:size(GCB_landflux_list,3)
    landflux_temp=GCB_landflux_list(:,:,i);
    GCB_landflux_anomaly(i,:)=GCB_landflux_2023_list(i,:)-mean(landflux_temp);
end
GCB_landflux_anomaly_std=std(GCB_landflux_anomaly);
Fire_anomaly=Fire_2023-nanmean(Fire_list);
GCB_NEP_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCB_landflux_anomaly_std,2));

% GPP距平均值
GPP_anomaly_mean=GPP_2023-nanmean(GPP_list);
for i=1:size(GPP_std_list,3)

    temp_list=GPP_std_list(:,:,i);
    temp2023=GPP_2023_list(i,:);
    GPP_anomaly_std(i,:)=temp2023-nanmean(temp_list);

end
for i=1:length(GPP_anomaly_std)
    GPP_temp=GPP_std_list(:,i);
    GPP_std(i)=std(bootstrp(1000,@mean,GPP_temp));
end
GPP_anomaly_std=GPP_std;

% ER距平均值与标准差
GCB_ER_anomaly_mean=GCB_ER_2023-nanmean(GCB_ER_list);
GCB_ER_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCB_landflux_anomaly_std,2)+power(GPP_anomaly_std,2));
%%
ff=figure
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
b=bar([-GCAS_landflux_2023,-GCB_satellite_landflux_2023,-GCB_obspack_landflux_2023;...
    -GCAS_Mean_fire_NEE_2023mean,-GCB_satellite_Mean_fire_NEE_2023mean,-GCB_obspack_Mean_fire_NEE_2023mean;...
    GCAS_ER_2023,GCB_satellite_ER_2023,GCB_obspack_ER_2023],'FaceColor','flat'); hold on
b(1).CData = [193,229,193]/255;
b(2).CData = [147,199,163]/255;
b(3).CData = [127,160,139]/255;

% %分组误差棒
[M,N]=size([-GCAS_landflux_2023,-GCB_satellite_landflux_2023,-GCB_obspack_landflux_2023;...
    -GCAS_Mean_fire_NEE_2023mean,-GCB_satellite_Mean_fire_NEE_2023mean,-GCB_obspack_Mean_fire_NEE_2023mean;...
    GCAS_ER_2023,GCB_satellite_ER_2023,GCB_obspack_ER_2023]);
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),[-GCAS_landflux_2023,-GCB_satellite_landflux_2023,-GCB_obspack_landflux_2023;...
    -GCAS_Mean_fire_NEE_2023mean,-GCB_satellite_Mean_fire_NEE_2023mean,-GCB_obspack_Mean_fire_NEE_2023mean;...
    GCAS_ER_2023,GCB_satellite_ER_2023,GCB_obspack_ER_2023],...
    [GCAS_landflux_2023std,GCB_satellite_landflux_2023std,GCB_obspack_landflux_2023std;...
    GCAS_Mean_fire_NEE_2023std,GCB_satellite_Mean_fire_NEE_2023std,GCB_obspack_Mean_fire_NEE_2023std;...
    GCAS_ER_2023std,GCB_satellite_ER_2023std,GCB_obspack_ER_2023std], ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

set(gca,'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'NBP','NEP','TER'},'FontName','Arial','fontsize',14)
% ylim([-0.8,0.8])
% set(gca,'YTick', [-0.8:0.4:0.8]);
text('string','a','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
lgd=legend({'GCAS-extra','GCB2024-satellite','GCB2024-\itin-situ'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','on','Location','best')

% lgd.ItemTokenSize = [15 8];
% lgd.Position=[0.632486079643154 0.0983566153034353 0.173876867962756 0.101283877365259]
ylabel('PgC','FontName','Arial','FontSize',14)
% title('Monthly flux anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')


nexttile
b=bar([-GCAS_landflux_2023anomaly_mean,-GCB_satellite_landflux_2023anomaly_mean,-GCB_obspack_landflux_2023anomaly_mean;...
    -GCAS_Mean_fire_NEE_2023anomaly_mean,-GCB_satellite_Mean_fire_NEE_2023anomaly_mean,-GCB_obspack_Mean_fire_NEE_2023anomaly_mean;...
    GCAS_ER_2023anomaly,GCB_satellite_ER_2023anomaly,GCB_obspack_ER_2023anomaly],'FaceColor','flat'); hold on
b(1).CData = [193,229,193]/255;
b(2).CData = [147,199,163]/255;
b(3).CData = [127,160,139]/255;

% %分组误差棒
[M,N]=size([-GCAS_landflux_2023anomaly_mean,-GCB_satellite_landflux_2023anomaly_mean,-GCB_obspack_landflux_2023anomaly_mean;...
    -GCAS_Mean_fire_NEE_2023anomaly_mean,-GCB_satellite_Mean_fire_NEE_2023anomaly_mean,-GCB_obspack_Mean_fire_NEE_2023anomaly_mean;...
    GCAS_ER_2023anomaly,GCB_satellite_ER_2023anomaly,GCB_obspack_ER_2023anomaly]);
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),[-GCAS_landflux_2023anomaly_mean,-GCB_satellite_landflux_2023anomaly_mean,-GCB_obspack_landflux_2023anomaly_mean;...
    -GCAS_Mean_fire_NEE_2023anomaly_mean,-GCB_satellite_Mean_fire_NEE_2023anomaly_mean,-GCB_obspack_Mean_fire_NEE_2023anomaly_mean;...
    GCAS_ER_2023anomaly,GCB_satellite_ER_2023anomaly,GCB_obspack_ER_2023anomaly],...
    [GCAS_landflux_2023anomaly_std,GCB_satellite_landflux_2023anomaly_std,GCB_obspack_landflux_2023anomaly_std;...
    GCAS_Mean_fire_NEE_2023anomaly_std,GCB_satellite_Mean_fire_NEE_2023anomaly_std,GCB_obspack_Mean_fire_NEE_2023anomaly_std;...
    GCAS_ER_2023anomaly_std,GCB_satellite_ER_2023anomaly_std,GCB_obspack_ER_2023anomaly_std], ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

set(gca,'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'NBP','NEP','TER'},'FontName','Arial','fontsize',14)
ylim([-0.8,0.8])
set(gca,'YTick', [-0.8:0.4:0.8]);
text('string','b','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
lgd=legend({'GCAS-extra','GCB2024-satellite','GCB2024-\itin-situ'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','on','Location','best')
% lgd.ItemTokenSize = [15 8];
% lgd.Position=[0.632486079643154 0.0983566153034353 0.173876867962756 0.101283877365259]
ylabel('PgC','FontName','Arial','FontSize',14)
% title('Monthly flux anomalies (GCASv2)','FontName','Arial','FontSize',14,'fontweight','bold')


n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
c = redblue();
colormap(n2,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (GCB2024-{\itin-situ})','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
colormap(n6,flipud(nclCM(399)));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('TER anomalies (GCB2024-{\itin-situ})','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

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

ylim([-0.17,0.23])
set(gca,'YTick', [-0.1:0.1:0.2]);
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','e','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
title('Monthly flux anomalies (GCB2024-{\itin-situ})','FontName','Arial','FontSize',14,'fontweight','bold')



set(gcf,'unit','pixels','position',[540,96,1202,1055]);
result=['E:\phd_file\Boreal_North_America\Result\V9\GCB2024_obspack_NEE_bar.png']
% print(result,ff,'-r600','-dpng' );