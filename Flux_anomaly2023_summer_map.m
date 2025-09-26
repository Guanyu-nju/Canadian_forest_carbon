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

for year=2015:2022

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\summer\NEE_summer_' num2str(year) '.tif']);
    NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    GPP_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_' num2str(year) '.tif']);
    GPP_mean(:,:,year-2014)=GPP_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\ER\total_carbon\summer\ER_summer_' num2str(year) '.tif']);
    ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
NEE_mean=nanmean(NEE_mean,3);
GPP_mean=nanmean(GPP_mean,3);
ER_mean=nanmean(ER_mean,3);



NEE2023=importdata("E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\summer\NEE_summer_2023.tif").*pixel_mask;
GPP2023=importdata("E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_2023.tif").*pixel_mask;
ER2023=importdata("E:\phd_file\Boreal_North_America\ER\total_carbon\summer\ER_summer_2023.tif").*pixel_mask;


NEP_2023Anomaly=-NEE2023+NEE_mean;
GPP_2023Anomaly=GPP2023-GPP_mean;
ER_2023Anomaly=ER2023-ER_mean;
%%
% % area ratio
% count_sum=sum(sum(~isnan(GCB_NEP_2023Anomaly)));
% count1=sum(sum(GCB_NEP_2023Anomaly>0));
% result=count1/(count_sum)
%% calculate mean value (GCB2024)

for year=2015:2022

    NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\summer\NEE_summer_' num2str(year) '.tif']);
    GCB_NEE_mean(:,:,year-2014)=NNE_mean_temp.*pixel_mask;

    ER_mean_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\summer\ER_summer_' num2str(year) '.tif']);
    GCB_ER_mean(:,:,year-2014)=ER_mean_temp.*pixel_mask;

end
GCB_NEE_mean=nanmean(GCB_NEE_mean,3);
GCB_ER_mean=nanmean(GCB_ER_mean,3);



GCB_NEE2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\summer\NEE_summer_2023.tif").*pixel_mask;
GCB_ER2023=importdata("E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\summer\ER_summer_2023.tif").*pixel_mask;


GCB_NEP_2023Anomaly=-GCB_NEE2023+GCB_NEE_mean;
GCB_ER_2023Anomaly=GCB_ER2023-GCB_ER_mean;
%%

count_sum=sum(sum(~isnan(GCB_NEP_2023Anomaly)));
count1=sum(sum(GCB_NEP_2023Anomaly>0));
result=count1/(count_sum)
%% GCAS
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;


for month=1:12


    for year=2015:2022

        % original fire
        GFAS_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFAS\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        GFED_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\GFED\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);
        Mean_fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year) '_' num2str(month) '.tif']);

        Fire_list(year-2014,month)=nansum(nansum(Mean_fire.*pixel_mask.*area_grid/(10^15)));
        NNE_mean_temp=importdata(['E:\phd_file\Boreal_North_America\Prior and posterior fluxes at 1deg resolution\Mean_value\total_carbon\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        NEE_list(year-2014,month)=nansum(nansum(NNE_mean_temp.*area_grid/(10^15)));

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
spring_NEP_anomaly_mean=nansum(NEP_anomaly_mean(3:5));
summer_NEP_anomaly_mean=nansum(NEP_anomaly_mean(6:8));
autumn_NEP_anomaly_mean=nansum(NEP_anomaly_mean(9:11));



Fire_anomaly=Fire_2023-nanmean(Fire_list);
spring_Fire_anomaly_mean=nansum(Fire_anomaly(3:5));
summer_Fire_anomaly_mean=nansum(Fire_anomaly(6:8));
autumn_Fire_anomaly_mean=nansum(Fire_anomaly(9:11));

for i=1:size(GCAS_landflux_list,3)
    landflux_temp=GCAS_landflux_list(:,:,i);
    GCAS_landflux_anomaly(i,:)=GCAS_landflux_2023_list(i,:)-mean(landflux_temp);
end
spring_GCAS_landflux_anomaly_std=std(nansum(GCAS_landflux_anomaly(:,3:5),2));
summer_GCAS_landflux_anomaly_std=std(nansum(GCAS_landflux_anomaly(:,6:8),2));
autumn_GCAS_landflux_anomaly_std=std(nansum(GCAS_landflux_anomaly(:,9:11),2));

spring_NEP_anomaly_std=sqrt(power(spring_Fire_anomaly_mean*0.2,2)+power(spring_GCAS_landflux_anomaly_std,2));
summer_NEP_anomaly_std=sqrt(power(summer_Fire_anomaly_mean*0.2,2)+power(summer_GCAS_landflux_anomaly_std,2));
autumn_NEP_anomaly_std=sqrt(power(autumn_Fire_anomaly_mean*0.2,2)+power(autumn_GCAS_landflux_anomaly_std,2));


% GPP距平均值
GPP_anomaly_mean=GPP_2023-nanmean(GPP_list);
spring_GPP_anomaly_mean=nansum(GPP_anomaly_mean(3:5));
summer_GPP_anomaly_mean=nansum(GPP_anomaly_mean(6:8));
autumn_GPP_anomaly_mean=nansum(GPP_anomaly_mean(9:11));

for i=1:size(GPP_std_list,3)

    temp_list=GPP_std_list(:,:,i);
    temp2023=GPP_2023_list(i,:);
    GPP_anomaly_std(i,:)=temp2023-nanmean(temp_list);

end

spring_GPP_anomaly_std=std(bootstrp(1000,@mean,sum(GPP_anomaly_std(:,3:5),2)));
summer_GPP_anomaly_std=std(bootstrp(1000,@mean,sum(GPP_anomaly_std(:,6:8),2)));
autumn_GPP_anomaly_std=std(bootstrp(1000,@mean,sum(GPP_anomaly_std(:,9:11),2)));

for i=1:length(GPP_anomaly_std)
    GPP_temp=GPP_std_list(:,i);
    GPP_std(i)=std(bootstrp(1000,@mean,GPP_temp));
end
GPP_anomaly_std=GPP_std;


ER_anomaly_mean=GCAS_ER_2023-nanmean(ER_list);
spring_ER_anomaly_mean=nansum(ER_anomaly_mean(3:5));
summer_ER_anomaly_mean=nansum(ER_anomaly_mean(6:8));
autumn_ER_anomaly_mean=nansum(ER_anomaly_mean(9:11));

spring_ER_std=sqrt(power(spring_Fire_anomaly_mean*0.2,2)+power(spring_GCAS_landflux_anomaly_std,2)+power(spring_GPP_anomaly_std,2));
summer_ER_std=sqrt(power(summer_Fire_anomaly_mean*0.2,2)+power(summer_GCAS_landflux_anomaly_std,2)+power(summer_GPP_anomaly_std,2));
autumn_ER_std=sqrt(power(autumn_Fire_anomaly_mean*0.2,2)+power(autumn_GCAS_landflux_anomaly_std,2)+power(autumn_GPP_anomaly_std,2));

anomaly_list=[spring_GPP_anomaly_mean,spring_ER_anomaly_mean,spring_NEP_anomaly_mean;summer_GPP_anomaly_mean,summer_ER_anomaly_mean,summer_NEP_anomaly_mean;...
    autumn_GPP_anomaly_mean,autumn_ER_anomaly_mean,autumn_NEP_anomaly_mean];
anomaly_std_list=[spring_GPP_anomaly_std,spring_ER_std,spring_NEP_anomaly_std;summer_GPP_anomaly_std,summer_ER_std,summer_NEP_anomaly_std;...
    autumn_GPP_anomaly_std,autumn_ER_std,autumn_NEP_anomaly_std];
%% GCB2024


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
spring_GCB_NEP_anomaly_mean=nansum(GCB_NEP_anomaly_mean(3:5));
summer_GCB_NEP_anomaly_mean=nansum(GCB_NEP_anomaly_mean(6:8));
autumn_GCB_NEP_anomaly_mean=nansum(GCB_NEP_anomaly_mean(9:11));

for i=1:size(GCB_landflux_list,3)
    landflux_temp=GCB_landflux_list(:,:,i);
    GCB_landflux_anomaly(i,:)=GCB_landflux_2023_list(i,:)-mean(landflux_temp);
end
spring_GCB_landflux_anomaly_std=std(nansum(GCB_landflux_anomaly(:,3:5),2));
summer_GCB_landflux_anomaly_std=std(nansum(GCB_landflux_anomaly(:,6:8),2));
autumn_GCB_landflux_anomaly_std=std(nansum(GCB_landflux_anomaly(:,9:11),2));

spring_GCB_NEP_anomaly_std=sqrt(power(spring_Fire_anomaly_mean*0.2,2)+power(spring_GCB_landflux_anomaly_std,2));
summer_GCB_NEP_anomaly_std=sqrt(power(summer_Fire_anomaly_mean*0.2,2)+power(summer_GCB_landflux_anomaly_std,2));
autumn_GCB_NEP_anomaly_std=sqrt(power(autumn_Fire_anomaly_mean*0.2,2)+power(autumn_GCB_landflux_anomaly_std,2));

% ER
GCB_ER_anomaly_mean=GCB_ER_2023-nanmean(GCB_ER_list);
spring_GCB_ER_anomaly_mean=nansum(GCB_ER_anomaly_mean(3:5));
summer_GCB_ER_anomaly_mean=nansum(GCB_ER_anomaly_mean(6:8));
autumn_GCB_ER_anomaly_mean=nansum(GCB_ER_anomaly_mean(9:11));

spring_GCB_ER_std=sqrt(power(spring_Fire_anomaly_mean*0.2,2)+power(spring_GCB_landflux_anomaly_std,2)+power(spring_GPP_anomaly_std,2));
summer_GCB_ER_std=sqrt(power(summer_Fire_anomaly_mean*0.2,2)+power(summer_GCB_landflux_anomaly_std,2)+power(summer_GPP_anomaly_std,2));
autumn_GCB_ER_std=sqrt(power(autumn_Fire_anomaly_mean*0.2,2)+power(autumn_GCB_landflux_anomaly_std,2)+power(autumn_GPP_anomaly_std,2));

GCB_anomaly_list=[spring_GPP_anomaly_mean,spring_GCB_ER_anomaly_mean,spring_GCB_NEP_anomaly_mean;summer_GPP_anomaly_mean,summer_GCB_ER_anomaly_mean,summer_GCB_NEP_anomaly_mean;...
    autumn_GPP_anomaly_mean,autumn_GCB_ER_anomaly_mean,autumn_GCB_NEP_anomaly_mean];
GCB_anomaly_std_list=[spring_GPP_anomaly_std,spring_GCB_ER_std,spring_GCB_NEP_anomaly_std;summer_GPP_anomaly_std,summer_GCB_ER_std,summer_GCB_NEP_anomaly_std;...
    autumn_GPP_anomaly_std,autumn_GCB_ER_std,autumn_GCB_NEP_anomaly_std];


Total_anomaly_list=[spring_GPP_anomaly_mean,spring_ER_anomaly_mean,spring_GCB_ER_anomaly_mean,spring_NEP_anomaly_mean,spring_GCB_NEP_anomaly_mean;...
    summer_GPP_anomaly_mean,summer_ER_anomaly_mean,summer_GCB_ER_anomaly_mean,summer_NEP_anomaly_mean,summer_GCB_NEP_anomaly_mean;...
    autumn_GPP_anomaly_mean,autumn_ER_anomaly_mean,autumn_GCB_ER_anomaly_mean,autumn_NEP_anomaly_mean,autumn_GCB_NEP_anomaly_mean];
Total_anomaly_std=[spring_GPP_anomaly_std,spring_ER_std,spring_GCB_ER_std,spring_NEP_anomaly_std,spring_GCB_NEP_anomaly_std;...
    summer_GPP_anomaly_std,summer_ER_std,summer_GCB_ER_std,summer_NEP_anomaly_std,summer_GCB_NEP_anomaly_std;...
    autumn_GPP_anomaly_std,autumn_ER_std,autumn_GCB_ER_std,autumn_NEP_anomaly_std,autumn_GCB_NEP_anomaly_std];
%% kNDVI
kNDVI2023=importdata("E:\phd_file\Boreal_North_America\kNDVI\summer\kNDVI_summer_2023.tif").*pixel_mask;
for year=2015:2022

    kNDVI_mean_temp=importdata(['E:\phd_file\Boreal_North_America\kNDVI\summer\kNDVI_summer_' num2str(year) '.tif']);
    kNDVI_mean(:,:,year-2014)=kNDVI_mean_temp.*pixel_mask;


end
kNDVI_mean=nanmean(kNDVI_mean,3);
kNDVI_2023Anomaly=kNDVI2023-kNDVI_mean;

for month=1:12


    for year=2015:2022
        kNDVI_mean_temp=importdata(['E:\phd_file\Boreal_North_America\kNDVI\month\globe\kNDVI_' num2str(year) '_' num2str(month) '.tif']);
        kNDVI_list(year-2014,month)=nansum(nansum(kNDVI_mean_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));
    end

end

for month=1:12
    year=2023;
    kNDVI_2023_temp=importdata(['E:\phd_file\Boreal_North_America\kNDVI\month\globe\kNDVI_' num2str(year) '_' num2str(month) '.tif']);
    kNDVI_2023(month)=nansum(nansum(kNDVI_2023_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));
end

kNDVI_anomaly_mean=kNDVI_2023-nanmean(kNDVI_list);
%% SIF
area_grid=importdata("E:\phd_file\Boreal_North_America\degree2meter.tif")*1000000.*pixel_mask;

SIF2023=importdata("E:\phd_file\Boreal_North_America\SIF\TROPOMI\summer\SIF_summer_2023.tif").*pixel_mask;
for year=2019:2022

    SIF_mean_temp=importdata(['E:\phd_file\Boreal_North_America\SIF\TROPOMI\summer\SIF_summer_' num2str(year) '.tif']);
    SIF_mean(:,:,year-2018)=SIF_mean_temp.*pixel_mask;


end
SIF_mean=nanmean(SIF_mean,3);
SIF_2023Anomaly=SIF2023-SIF_mean;

for month=1:12


    for year=2019:2022
        SIF_mean_temp=importdata(['E:\phd_file\Boreal_North_America\SIF\TROPOMI\month\1degree\SIF_' num2str(year) '_' num2str(month) '.tif']);
        SIF_list(year-2018,month)=nansum(nansum(SIF_mean_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));
    end

end

for month=1:12
    year=2023;
    SIF_2023_temp=importdata(['E:\phd_file\Boreal_North_America\SIF\TROPOMI\month\1degree\SIF_' num2str(year) '_' num2str(month) '.tif']);
    SIF_2023(month)=nansum(nansum(SIF_2023_temp.*pixel_mask.*area_grid))/(nansum(nansum(area_grid)));
end

SIF_anomaly_mean=SIF_2023-nanmean(SIF_list);
%%
% set(0, 'DefaultFigureRenderer', 'painters');
f=figure
t = tiledlayout(3,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';


% part 1 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n2,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} summer^{-1}','fontsize',14);
title('NEP anomalies (GCAS-extra)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','a','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 2 *********************************************************************************
n2=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_NEP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
c = redblue();
colormap(n2,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} summer^{-1}','fontsize',14);
title('NEP anomalies (GCB2024-satellite)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','b','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')


% part 3 *********************************************************************************
n4=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GPP_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n4,nclCM(399));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} summer^{-1}','fontsize',14);
title('GPP anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','c','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 4 *********************************************************************************
n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n6,flipud(nclCM(399)));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} summer^{-1}','fontsize',14);
title('TER anomalies (GCAS-extra)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','d','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 5 *********************************************************************************
n6=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,GCB_ER_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n6,flipud(nclCM(399)));

caxis([-120,120]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-120:60:120]);
set(get(h,'ylabel'),'string','gC m^{-2} summer^{-1}','fontsize',14);
title('TER anomalies (GCB2024-satellite)','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','e','Units','normalized','position',[-0.155950217182262 1.15711095257911 0],'FontName','Arial','FontSize',20,'fontweight','bold')

% part 6 *********************************************************************************
n7=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,kNDVI_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n7,nclCM(399));

caxis([-0.06,0.06]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.06:0.03:0.06]);
set(get(h,'ylabel'),'string','kNDVI','fontsize',14);
title('kNDVI anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','f','Units','normalized','position',[-0.160574295917181 1.08988406182281 0],'FontName','Arial','FontSize',20,'fontweight','bold')


% part 7 ************************************************************************<0*********
nexttile
b=bar(Total_anomaly_list,'FaceColor','flat'); hold on
b(1).CData = [17,119,51]/255;
b(2).CData = [141,124,106]/255;
b(3).CData = [200,184,147]/255;
b(4).CData = [120,149,203]/255;
b(5).CData = [184,128,197]/255;


[M,N]=size(Total_anomaly_list);
for i=1:N
    xx(:,i)=b(i).XEndPoints';
end

h2=errorbar(xx(:,:),Total_anomaly_list(:,:),Total_anomaly_std(:,:), ...
    'LineStyle', 'none', 'Color', 'k', 'LineWidth',0.8,'CapSize',7,'Marker','none')

ylim([-0.5,0.5])
set(gca,'YTick', [-0.4:0.2:0.4]);
ylabel('PgC season^{-1}','FontName','Arial','FontSize',14);
xlabel('Season','FontName','Arial','FontSize',14)
set(gca,'XTick',[1:1:3],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'MAM','JJA','SON'},'FontName','Arial','fontsize',14)
text('string','g','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
lgd=legend({'GPP','TER (GCAS-extra)','TER (GCB2024-satellite)','NEP (GCAS-extra)','NEP (GCB2024-satellite)'},'NumColumns',1,'FontName','Arial','FontSize',10,'Box','on','Location','southeast')
lgd.ItemTokenSize = [9 9];
lgd.Position=[0.160277139365076 0.0752247395767492 0.0987297896011873 0.0927723814625828]
% title('Seasonal flux anomalies','FontName','Arial','FontSize',14,'fontweight','bold')

% part 8 ************************************************************************<0*********
nexttile
yyaxis left
plot([0.5,9.5],[0,0],'--','LineWidth',1,'MarkerSize',3,'color',[242,123,119]/255);hold on
p1=plot(SIF_anomaly_mean(3:11),'-','LineWidth',3,'MarkerSize',15,'color',[242,123,119]/255);
ylim([-0.035,0.06])
set(gca,'YTick', [-0.03:0.03:0.06]);
ylabel('SIF (mW m^{-2} nm^{-1} sr^{-1})','FontName','Arial','FontSize',14);
set(gca,'YColor',[242,123,119]/255)

yyaxis right
p2=plot(kNDVI_anomaly_mean(3:11),'-','LineWidth',3,'MarkerSize',15,'color',[246,173,60]/255);hold on
ylim([-0.035,0.06])
set(gca,'YTick', [-0.03:0.03:0.06]);
ylabel('kNDVI','FontName','Arial','FontSize',14);
set(gca,'YColor',[246,173,60]/255)

% ylabel('PgC month^{-1}','FontName','Arial','FontSize',14);
xlim([0.5,9.5])
set(gca,'XTick',[1:1:9],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'M','A','M','J','J','A','S','O','N'},'FontName','Arial','fontsize',14)
text('string','h','Units','normalized','position',[-0.158500470010452 1.08805282109315 0],'FontName','Arial','FontSize',20,'fontweight','bold')
legend([p1 p2],{'SIF','kNDVI'},'NumColumns',1,'FontName','Arial','FontSize',14,'Box','off','Location','northeast')
xlabel('Month','FontName','Arial','FontSize',14)
% title('Monthly flux anomalies','FontName','Arial','FontSize',14,'fontweight','bold')

% part 9 *********************************************************************************
n9=nexttile
m_proj('lambert','long',[-150 -50],'lat',[41 75]);
m_patch(bou_canX,bou_canY,[190 190 190]/255,'linewidth',0.5); hold on
% m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',0.7); hold on
m_pcolor(lon1,lat1,SIF_2023Anomaly,'linestyle','none');
% m_plot(bou_canX,bou_canY,'linewidth',0.1,'color','k');
m_grid('tickdir','in','FontName','Arial','xtick',[-150:30:-40],'ytick',[45:10:65],...
    'linewidth',0.7,'xaxisloc','bottom','yaxisloc','left','FontSize',14,'xticklabels',[-160:30:-40],'yticklabels',[45:10:65]);
colormap(n9,nclCM(399));

caxis([-0.08,0.08]);
h=colorbar('location','eastoutside','FontName','Arial','FontSize',14,'ytick',[-0.08:0.04:0.08]);
set(get(h,'ylabel'),'string','mW m^{-2} sr^{-1} nm^{-1}','fontsize',14);
title('SIF anomalies','FontName','Arial','FontSize',14,'fontweight','bold')
text('string','i','Units','normalized','position',[-0.160574295917181 1.08988406182281 0],'FontName','Arial','FontSize',20,'fontweight','bold')

set(gcf,'unit','centimeters','position',[10.054166666666667,8.128,45.825833333333335,24.510999999999996]);
result=['E:\phd_file\Boreal_North_America\Result\V9\Flux_anomaly2023_summer_map.png']
% print(result,f,'-r600','-dpng');

