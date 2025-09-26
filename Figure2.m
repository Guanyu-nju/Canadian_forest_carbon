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
%% calculate mean value
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;

for year=2015:2022

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\year\NEE_' num2str(year) '.tif']);
    NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\year\ER_' num2str(year) '.tif']);
    ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
NEE_mean=nanmean(NEE_mean,3);
GPP_mean=nanmean(GPP_mean,3);
ER_mean=nanmean(ER_mean,3);



NEE2023=importdata("E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\year\NEE_2023.tif").*pixel_mask;
GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_2023.tif").*pixel_mask;
ER2023=importdata("E:\phd_file\Boreal_North_America\ER\total_carbon\year\ER_2023.tif").*pixel_mask;


NEP_2023Anomaly=-NEE2023+NEE_mean;
GPP_2023Anomaly=GPP2023-GPP_mean;
ER_2023Anomaly=ER2023-ER_mean;
%%
% 面积占比
count_sum=sum(sum(~isnan(NEP_2023Anomaly)));
count1=sum(sum(NEP_2023Anomaly>0));
result=count1/(count_sum)
%% calculate mean value (GCB2024)

for year=2015:2022

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\year\NEE_' num2str(year) '.tif']);
    GCB_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);
    GCB_GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\year\ER_' num2str(year) '.tif']);
    GCB_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);
GCB_GPP_mean=nanmean(GCB_GPP_mean,3);
GCB_ER_mean=nanmean(GCB_ER_mean,3);



GCB_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\year\NEE_2023.tif").*pixel_mask;
GCB_GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_2023.tif").*pixel_mask;
GCB_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\year\ER_2023.tif").*pixel_mask;


GCB_NEP_2023Anomaly=-GCB_NEE2023+GCB_NEE_mean;
GPP_2023Anomaly=GCB_GPP2023-GCB_GPP_mean;
GCB_ER_2023Anomaly=GCB_ER2023-GCB_ER_mean;
%%
% area ratio
count_sum=sum(sum(~isnan(GCB_NEP_2023Anomaly)));
count1=sum(sum(GCB_NEP_2023Anomaly>0));
result=count1/(count_sum)
%% GCAS

% time list
for month=1:12


    for year=2015:2022

        % original fire
        GFAS_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        GFED_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Fire_list(year-2014,month)=nansum(nansum(Mean_fire.*pixel_mask.*area_grid/(10^15)));

        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));

        BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_total_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_total_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_total_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
        CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_total_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);

        % original landflux
        BEPS_GFAS_landflux=BEPS_GFAS_NEE+GFAS_fire;
        BEPS_GFED_landflux=BEPS_GFED_NEE+GFED_fire;
        CASA_GFAS_landflux=CASA_GFAS_NEE+GFAS_fire;
        CASA_GFED_landflux=CASA_GFED_NEE+GFED_fire;
        GCAS_landflux_list(year-2014,month,1)=nansum(nansum(BEPS_GFAS_landflux.*area_grid/(10^15)));
        GCAS_landflux_list(year-2014,month,2)=nansum(nansum(BEPS_GFED_landflux.*area_grid/(10^15)));
        GCAS_landflux_list(year-2014,month,3)=nansum(nansum(CASA_GFAS_landflux.*area_grid/(10^15)));
        GCAS_landflux_list(year-2014,month,4)=nansum(nansum(CASA_GFED_landflux.*area_grid/(10^15)));


        GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
        GPP_list(year-2014,month)=nansum(nansum(GPP_mean_temp.*pixel_mask.*area_grid/(10^15)));

        GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
        Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

        GPP_std_list(year-2014,month,1)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
        GPP_std_list(year-2014,month,2)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));


        ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\month\ER_' num2str(year) '_' num2str(month) '.tif']);
        ER_list(year-2014,month)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));


    end

end

% 2023 data list
for month=1:12


    year=2023;
    % original fire
    GFAS_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        GFED_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
    Fire_2023(month)=nansum(nansum(Mean_fire.*pixel_mask.*area_grid/(10^15)));

    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    NEE_2023(month)=nansum(nansum(NNE_2023_temp.*pixel_mask.*area_grid/(10^15)));

    BEPS_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFAS\opt_total_monthly\Opt_NEE_BEPS_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    BEPS_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.BEPS_GFED\opt_total_monthly\Opt_NEE_BEPS_GFED_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFAS_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFAS\opt_total_monthly\Opt_NEE_CASA_GFAS_' num2str(year) '_' num2str(month) '.tif']);
    CASA_GFED_NEE=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\posterior.fluxes.CASA_GFED\opt_total_monthly\Opt_NEE_CASA_GFED_' num2str(year) '_' num2str(month) '.tif']);

    % original landflux
    BEPS_GFAS_landflux=BEPS_GFAS_NEE+GFAS_fire;
    BEPS_GFED_landflux=BEPS_GFED_NEE+GFED_fire;
    CASA_GFAS_landflux=CASA_GFAS_NEE+GFAS_fire;
    CASA_GFED_landflux=CASA_GFED_NEE+GFED_fire;

    GCAS_landflux_2023_list(1,month)=nansum(nansum(BEPS_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_2023_list(2,month)=nansum(nansum(BEPS_GFED_landflux.*area_grid/(10^15)));
    GCAS_landflux_2023_list(3,month)=nansum(nansum(CASA_GFAS_landflux.*area_grid/(10^15)));
    GCAS_landflux_2023_list(4,month)=nansum(nansum(CASA_GFED_landflux.*area_grid/(10^15)));


    GPP_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);
    GPP_2023(month)=nansum(nansum(GPP_2023_temp.*pixel_mask.*area_grid/(10^15)));

    GOSIF_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\GOSIF\monthly\GPP_GOSIF_' num2str(year) '_' num2str(month) '.tif']);
    Fluxsat_GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\FluxSat\monthly\GPP_FluxSat_' num2str(year) '_' num2str(month) '.tif']);

    GPP_2023_list(1,month)=nansum(nansum(GOSIF_GPP_temp.*pixel_mask.*area_grid/(10^15)));
    GPP_2023_list(2,month)=nansum(nansum(Fluxsat_GPP_temp.*pixel_mask.*area_grid/(10^15)));


    ER_2023_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\month\ER_' num2str(year) '_' num2str(month) '.tif']);
    GCAS_ER_2023(month)=nansum(nansum(ER_2023_temp.*pixel_mask.*area_grid/(10^15)));


end

% 2023 flux anomaly
NEP_anomaly_mean=-NEE_2023+nanmean(NEE_list);
Fire_anomaly=Fire_2023-nanmean(Fire_list);

for i=1:size(GCAS_landflux_list,3)
    landflux_temp=GCAS_landflux_list(:,:,i);
    GCAS_landflux_anomaly(i,:)=GCAS_landflux_2023_list(i,:)-mean(landflux_temp);
end
GCAS_landflux_anomaly_std=std(GCAS_landflux_anomaly);
NEP_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCAS_landflux_anomaly_std,2));

% GPP anomaly
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

% ERanomaly and std
ER_anomaly_mean=GCAS_ER_2023-nanmean(ER_list);

ER_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCAS_landflux_anomaly_std,2)+power(GPP_anomaly_std,2));

%% GCB2024

% time series
for month=1:12


    for year=2015:2022

        Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);

        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GCB_NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*pixel_mask.*area_grid/(10^15)));


        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year) '_' num2str(month) '.tif']);

        % original landflux
        A_landflux=A_NEE+Mean_fire;
        B_landflux=B_NEE+Mean_fire;
        C_landflux=C_NEE+Mean_fire;
        D_landflux=D_NEE+Mean_fire;
        E_landflux=E_NEE+Mean_fire;
        F_landflux=F_NEE+Mean_fire;

        GCB_landflux_list(year-2014,month,1)=nansum(nansum(A_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,2)=nansum(nansum(B_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,3)=nansum(nansum(C_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,4)=nansum(nansum(D_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,5)=nansum(nansum(E_landflux.*area_grid/(10^15)));
        GCB_landflux_list(year-2014,month,6)=nansum(nansum(F_landflux.*area_grid/(10^15)));

        ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']);
        GCB_ER_list(year-2014,month)=nansum(nansum(ER_mean_temp.*pixel_mask.*area_grid/(10^15)));

    end

end

% 2023 year data list
for month=1:12


    year=2023;
    Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
    NNE_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
    GCB_NEE_2023(month)=nansum(nansum(NNE_2023_temp.*pixel_mask.*area_grid/(10^15)));


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year) '_' num2str(month) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year) '_' num2str(month) '.tif']);

    % original landflux
    A_landflux=A_NEE+Mean_fire;
    B_landflux=B_NEE+Mean_fire;
    C_landflux=C_NEE+Mean_fire;
    D_landflux=D_NEE+Mean_fire;
    E_landflux=E_NEE+Mean_fire;
    F_landflux=F_NEE+Mean_fire;

    GCB_landflux_2023_list(1,month)=nansum(nansum(A_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(2,month)=nansum(nansum(B_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(3,month)=nansum(nansum(C_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(4,month)=nansum(nansum(D_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(5,month)=nansum(nansum(E_landflux.*pixel_mask.*area_grid/(10^15)));
    GCB_landflux_2023_list(6,month)=nansum(nansum(F_landflux.*pixel_mask.*area_grid/(10^15)));

    ER_2023_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']);
    GCB_ER_2023(month)=nansum(nansum(ER_2023_temp.*pixel_mask.*area_grid/(10^15)));


end

% 2023 flux anomaly
GCB_NEP_anomaly_mean=-GCB_NEE_2023+nanmean(GCB_NEE_list);

for i=1:size(GCB_landflux_list,3)
    landflux_temp=GCB_landflux_list(:,:,i);
    GCB_landflux_anomaly(i,:)=GCB_landflux_2023_list(i,:)-mean(landflux_temp);
end
GCB_landflux_anomaly_std=std(GCB_landflux_anomaly);
GCB_NEP_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCB_landflux_anomaly_std,2));

% ER anomaly and encertainty
GCB_ER_anomaly_mean=GCB_ER_2023-nanmean(GCB_ER_list);
GCB_ER_anomaly_std=sqrt(power(Fire_anomaly*0.2,2)+power(GCB_landflux_anomaly_std,2)+power(GPP_anomaly_std,2));


%%
for year=2015:2023

    Con_April=importdata(['E:\phd_file\Boreal_North_America\OCO2-XCO2\OCO2-XCO2_' num2str(year) '_4.tif']);
    April_mask=Con_April;
    April_mask(~isnan(April_mask))=1;
    Con_August=importdata(['E:\phd_file\Boreal_North_America\OCO2-XCO2\OCO2-XCO2_' num2str(year) '_8.tif']);

    April_mask(isnan(Con_August))=nan;
    Con_diff=nansum(nansum(Con_April.*pixel_mask.*area_grid.*April_mask))/(nansum(nansum(area_grid.*April_mask)))-nansum(nansum(Con_August.*pixel_mask.*area_grid.*April_mask))/nansum(nansum(area_grid.*April_mask));


    Data_list(year-2014)=Con_diff;


end
%% Atmospheric site
etl_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_etl_monthly_co2_SG.xlsx");
result = innerjoin(etl_data(etl_data.months == 4, :), etl_data(etl_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
etl_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

inu_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_inu_monthly_co2_SG.xlsx");
result = innerjoin(inu_data(inu_data.months == 4, :), inu_data(inu_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
inu_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

fsd_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_fsd_monthly_co2_SG.xlsx");
result = innerjoin(fsd_data(fsd_data.months == 4, :), fsd_data(fsd_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
fsd_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

est_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_est_monthly_co2_SG.xlsx");
result = innerjoin(est_data(est_data.months == 4, :), est_data(est_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
est_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

egb_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_egb_monthly_co2_SG.xlsx");
result = innerjoin(egb_data(egb_data.months == 4, :), egb_data(egb_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
egb_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

bra_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_bra_monthly_co2_SG.xlsx");
result = innerjoin(bra_data(bra_data.months == 4, :), bra_data(bra_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
bra_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

bck_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_bck_monthly_co2_SG.xlsx");
result = innerjoin(bck_data(bck_data.months == 4, :), bck_data(bck_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
bck_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

wsa_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_wsa_monthly_co2_SG.xlsx");
result = innerjoin(wsa_data(wsa_data.months == 4, :), wsa_data(wsa_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
wsa_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

dws_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\Daytime_dws_monthly_co2_SG.xlsx");
result = innerjoin(dws_data(dws_data.months == 4, :), dws_data(dws_data.months == 8, :), 'Keys', 'years');
co2_diff=result{:,4}-result{:,7};
dws_co2_diff=co2_diff(end)-nanmean(co2_diff(1:end-1)); clear result co2_diff

site_co2_diff=[etl_co2_diff,inu_co2_diff,fsd_co2_diff,egb_co2_diff,bra_co2_diff,bck_co2_diff,wsa_co2_diff,dws_co2_diff];
lon_co2_diff=[-104.9868,-133.5342,-81.57,-79.7838,-104.7113,-115.9194,-60.0093,-79.468];
lat_co2_diff=[54.3541,68.3178,49.8752,44.231,50.2017,62.7981,43.9322,43.7805];

%% EC site
% SOB
SOB_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\SOB_monthly.xlsx"); clear year
SOB_data.Year = year(SOB_data.Time_UTC);
SOB_annual_sums = varfun(@sum, SOB_data, 'GroupingVariables', 'Year', ...
                     'InputVariables', {'EcosystemRespiration','GrossEcosystemPhotosynthesis', 'NetEcosystemProduction'});
SOB_NEP_2023anomaly=SOB_annual_sums{end,5}-nanmean(SOB_annual_sums{1:end-1,5});
SOB_GPP_2023anomaly=SOB_annual_sums{end,4}-nanmean(SOB_annual_sums{1:end-1,4});
SOB_ER_2023anomaly=SOB_annual_sums{end,3}-nanmean(SOB_annual_sums{1:end-1,3});

% OJP
OJP_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\OJP_monthly.xlsx"); clear year
OJP_data.Year = year(OJP_data.Time_UTC);
OJP_annual_sums = varfun(@sum, OJP_data, 'GroupingVariables', 'Year', ...
                     'InputVariables', {'EcosystemRespiration','GrossEcosystemPhotosynthesis', 'NetEcosystemProduction'});
OJP_NEP_2023anomaly=OJP_annual_sums{end,5}-nanmean(OJP_annual_sums{1:end-1,5});
OJP_GPP_2023anomaly=OJP_annual_sums{end,4}-nanmean(OJP_annual_sums{1:end-1,4});
OJP_ER_2023anomaly=OJP_annual_sums{end,3}-nanmean(OJP_annual_sums{1:end-1,3});

% SCC
SCC_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\Guanyu_ForestSites_2015-2023.xlsx",'Sheet', 'SCC');
SCC_annual_sums = varfun(@sum, SCC_data, 'GroupingVariables', 'year', ...
                     'InputVariables', {'nee','gpp', 'reco'});
SCC_annual_sums{6,5}=nan;
SCC_annual_sums{:,3}=-SCC_annual_sums{:,3};
SCC_annual_sums{:,4}=-SCC_annual_sums{:,4};
SCC_NEP_2023anomaly=SCC_annual_sums{end,3}-nanmean(SCC_annual_sums{1:end-1,3});
SCC_GPP_2023anomaly=SCC_annual_sums{end,4}-nanmean(SCC_annual_sums{1:end-1,4});
SCC_ER_2023anomaly=SCC_annual_sums{end,5}-nanmean(SCC_annual_sums{1:end-1,5});

% HPC
HPC_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\Guanyu_ForestSites_2015-2023.xlsx",'Sheet', 'HPC');
HPC_annual_sums = varfun(@sum, HPC_data, 'GroupingVariables', 'year', ...
                     'InputVariables', {'nee','gpp', 'reco'});

HPC_annual_sums{:,3}=-HPC_annual_sums{:,3};
HPC_annual_sums{:,4}=-HPC_annual_sums{:,4};
HPC_NEP_2023anomaly=HPC_annual_sums{end,3}-nanmean(HPC_annual_sums{1:end-1,3});
HPC_GPP_2023anomaly=HPC_annual_sums{end,4}-nanmean(HPC_annual_sums{1:end-1,4});
HPC_ER_2023anomaly=HPC_annual_sums{end,5}-nanmean(HPC_annual_sums{1:end-1,5});

% TPD
for year=2015:2023
    data_temp=readtable(['E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\CA-TPD\TPD_Monthly_' num2str(year) '.csv']);
    TPD_annual_sum(year-2014,1)=sum(data_temp{1:12,3});
    TPD_annual_sum(year-2014,2)=sum(data_temp{1:12,4});
    TPD_annual_sum(year-2014,3)=sum(data_temp{1:12,5});
end
TPD_NEP_2023anomaly=TPD_annual_sum(end,3)-mean(TPD_annual_sum(1:8,3));

% TP3
for year=2015:2023
    data_temp=readtable(['E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\CA-TPD\CA_TP3_Monthly_' num2str(year) '.csv']);
    TP3_annual_sum(year-2014,1)=sum(data_temp{1:12,4})+sum(data_temp{1:12,5});
    TP3_annual_sum(year-2014,2)=sum(data_temp{1:12,4});
    TP3_annual_sum(year-2014,3)=sum(data_temp{1:12,5});

end
TP3_NEP_2023anomaly=TP3_annual_sum(end,3)-mean(TP3_annual_sum(1:8,3))

EC_NEP_2023anomaly=[SOB_NEP_2023anomaly,SCC_NEP_2023anomaly,HPC_NEP_2023anomaly,TP3_NEP_2023anomaly,TPD_NEP_2023anomaly];
lon_EC=[-105.118,-121.2990,-133.5188,-80.3483,-80.5577];
lat_EC=[53.987,61.3082,68.3203,42.7068,42.6353];


%%
% set(0, 'DefaultFigureRenderer', 'painters');
f=figure
t = tiledlayout(3,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','centimeters','position',[10.054166666666669,8.128,46.646041666666676,24.510999999999996]);


% part 1 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% shapefile
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% shapefile
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %grid 
c = redblue();
colormap(n2,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (GCAS-extra)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 2 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% shapefile
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% shapefile
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %grid
c = redblue();
colormap(n2,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (GCB2024-satellite)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')


% part 3 *********************************************************************************
n3=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on % 填充矢量边界
m_scatter(lon_EC,lat_EC,80,EC_NEP_2023anomaly','filled','o','MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k','linewi',0.5);
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
colormap(n3,nclCM(399));
caxis([-40,40]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('NEP anomalies (EC observations)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
% text('string','CA-SCC','Units','normalized','position',[0.383166228587275 0.565583156518392 0],'FontName','Arial','FontSize',10,'fontweight','bold')


% part 4 *********************************************************************************
n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
colormap(n6,flipud(nclCM(399)));
caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('TER anomalies (GCAS-extra)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 5 *********************************************************************************
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
title('TER anomalies (GCB2024-satellite)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','e','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')


% part 6 *********************************************************************************
n4=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on% 填充矢量边界
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GPP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');% 填充矢量边界
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
colormap(n4,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} yr^{-1}','fontsize',14);
title('GPP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','f','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')


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


ylim([-0.17,0.23])
set(gca,'YTick', [-0.1:0.1:0.2]);
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','g','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
title('Monthly flux anomalies (GCAS-extra)','FontName','Arial','FontSize',14,'fontweight','bold')

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

ylim([-0.17,0.23])
set(gca,'YTick', [-0.1:0.1:0.2]);
ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,12.5])
set(gca,'XTick',[1:1:12],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',},'FontName','Arial','fontsize',14)
text('string','h','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
legend([p1 p2,p3],{'NEP','GPP','TER'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
title('Monthly flux anomalies (GCB2024-satellite)','FontName','Arial','FontSize',14,'fontweight','bold')

n7=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on % 填充矢量边界
m6=m_scatter(lon_co2_diff,lat_co2_diff,80,site_co2_diff,'filled','^','MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k','linewi',0.5);
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
colormap(n7,flipud(nclCM(148)));
caxis([-4,4]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-4:2:4]);
set(get(h,'ylabel'),'string','ppm','fontsize',14);
title('Seasonal CO_2 amplitude anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','i','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')

n10=axes('position',[0.672490074314944 0.0699481865284974 0.083607486660666 0.105039660603392])
b1=bar(Data_list,'FaceColor','flat');hold on
for i=1:8
    b1.CData(i,:) = [55,140,230]/255;
end
b1.CData(9,:) = [237,145,29]/255;
ylim([6,10.5])
set(gca,'YTick', [6:2:10]);
ylabel('ppm','FontName','Arial','FontSize',10);
xlim([0.25,9.75])
set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',10)
set(gca,'XTickLabel',[2015:1:2023],'FontName','Arial','fontsize',10)
xtickangle(90)
Data_list(end)-nanmean(Data_list(1:end-1))
% n11=axes('position',[0.672490074314944 0.700142233059027 0.083607486660666 0.105039660603392])
% b2=bar([SCC_NEP_2023anomaly,SCC_GPP_2023anomaly,SCC_ER_2023anomaly],'FaceColor',[211,213,212]/255);hold on
% ylim([-90,35])
% set(gca,'YTick', [-90:30:30]);
% % ylabel('Flux anomalies','FontName','Arial','FontSize',10);
% % xlim([0.25,9.75])
% % set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',12)
% set(gca,'XTickLabel',{'NEP','GPP','TER'},'FontName','Arial','fontsize',10)
% text('string','CA-SCC','Units','normalized','position',[0.5897988816485 0.871557975732325 0],'FontName','Arial','FontSize',10,'fontweight','bold')



% n10=axes('position',[0.709926262630319,0.395704867210711,0.222310456184476,0.232650855568847],'color','none')
% m_proj('lambert','long',[-150 -50],'lat',[41 75]);
% % m_scatter(lon_EC,lat_EC,40,EC_NEP_2023anomaly,'filled','^','MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k','linewi',0.5);
% % m6=m_scatter(lon_co2_diff,lat_co2_diff,120,site_co2_diff','filled','^','MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k','linewi',0.5);
% 
% m_scatter(lon_EC,lat_EC,40,EC_NEP_2023anomaly','filled','o','MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'k','linewi',0.5);
% m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
%     'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]); %对齐格网
% % colormap(n10,hot);
% 
% colormap(n10,nclCM(399,8));
% caxis([-40,40]);
% set(n10,'clipping','off','layer','top'),axis off, hidden off,box off
% 
% h2 = colorbar(n10,'FontName','Arial','FontSize',14,'ytick',[-40:20:40]);
% set(get(h2,'title'),'string','gC m^{-2} yr^{-1}','fontsize',14);
% set(h2,'position',[0.936018150879182 0.392314335060449 0.00907543959160495 0.232650855568847])


set(gcf,'unit','centimeters','position',[10.054166666666669,8.128,46.646041666666676,24.510999999999996]);
result=['E:\phd_file\Boreal_North_America\Result\V9\Flux_Aomalies_year_map.png']
% print(result,f,'-r600','-dpng');

