%% CarboScope month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,1);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CarboScope (Prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,1);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\Prior\month\CarboScope_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CarboScope year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\year\CarboScope_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CAMS (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,2);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS\Prior\month\CAMS_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CAMS month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,2);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CAMS year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS\year\CAMS_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CTE (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,3);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CTE\Prior\month\CTE_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CTE month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,3);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CTE year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CTE\year\CTE_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% NISMON (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,4);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\prior\month\NISMON_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% NISMON month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,4);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% NISMON year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\year\NISMON_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CT (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,5);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\prior\month\CT_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CT month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,5);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CT year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\year\CT_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CMSF (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,6);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\prior\month\CMSF_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CMSF month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,6);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CMSF year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CMSF summer
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\summer\CMSF_NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CAMSS (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,7);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\prior\month\CAMSS_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CAMSS month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,7);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% CAMSS year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% CAMSS summer
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\summer\CAMSS_NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% GONGGA (Prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,8);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\prior\month\GONGGA_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% GONGGA month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,8);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% GONGGA year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% GONGGA summer
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\summer\GONGGA_NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% COLA month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,9);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% COLA year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\COLA\year\COLA_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% GCASv2 (Prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,10);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\Prior\month\GCASv2_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% GCASv2 month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,10);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% GCASv2 year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% GCASv2 summer
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\summer\GCASv2_NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% UoE month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,11);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% UoE year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\UoE\year\UoE_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% IAPCAS month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,12);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% IAPCAS year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\year\IAPCAS_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% MIROC-ACTM (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,13);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\prior\month\MIROC-ACTM_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% MIROC-ACTM month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,13);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% MIROC-ACTM year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\year\MIROC-ACTM_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% NTFVAR (prior) month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'prior_flux_land');
Boreal_data=Boreal_data(:,:,25:132,14);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\prior\month\NTFVAR_prior_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% NTFVAR month
clc
clear
Boreal_data=ncread("E:\phd_file\Boreal_North_America\GCB2024\GCP2024_inversions_1x1_version1_1_20240912.nc",'land_flux_only_fossil_cement_adjusted');
Boreal_data=Boreal_data(:,:,25:132,14);

for year=1:9
    for month=1:12

        fire=importdata(['E:\phd_file\Boreal_North_America\fire emission\mean\total_carbon\month\Fire_' num2str(year+2014) '_' num2str(month) '.tif']);
        Boreal_temp = Boreal_data(:,:,(year-1)*12+month);
        Boreal_temp=double(Boreal_temp)*10^15/12;
        Boreal_temp(Boreal_temp==0)=nan;
        Boreal_temp=rot90(Boreal_temp);
        Boreal_temp=Boreal_temp-fire;
        filename=['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year+2014) '_' num2str(month) '.tif']
        wtif(filename,Boreal_temp,[-180,180],[-90,90] );

    end
end
%% NTFVAR year
clc
clear

for year=2015:2023

    for month=1:12
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\year\NTFVAR_NEE_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% mean year

for year=2015:2023


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\year\CAMS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\year\CarboScope_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\year\COLA_NEE_' num2str(year) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\year\CTE_NEE_' num2str(year) '.tif']);
    G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\year\CT_NEE_' num2str(year) '.tif']);
    H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']);
    I_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']);
    J_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\year\IAPCAS_NEE_' num2str(year) '.tif']);
    K_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\year\MIROC-ACTM_NEE_' num2str(year) '.tif']);
    L_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\year\NISMON_NEE_' num2str(year) '.tif']);
    M_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\year\NTFVAR_NEE_' num2str(year) '.tif']);
    N_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\year\UoE_NEE_' num2str(year) '.tif']);



    NEE_mean(:,:,1)=A_NEE;
    NEE_mean(:,:,2)=B_NEE;
    NEE_mean(:,:,3)=C_NEE;
    NEE_mean(:,:,4)=D_NEE;
    NEE_mean(:,:,5)=E_NEE;
    NEE_mean(:,:,6)=F_NEE;
    NEE_mean(:,:,7)=G_NEE;
    NEE_mean(:,:,8)=H_NEE;
    NEE_mean(:,:,9)=I_NEE;
    NEE_mean(:,:,10)=J_NEE;
    NEE_mean(:,:,11)=K_NEE;
    NEE_mean(:,:,12)=L_NEE;
    NEE_mean(:,:,13)=M_NEE;
    NEE_mean(:,:,14)=N_NEE;

    NEE_median=nanmedian(NEE_mean,3);
    NEE_mean=nanmean(NEE_mean,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_mean,[-180,180],[-90,90] );
    filename=['E:\phd_file\Boreal_North_America\GCB2024\Median_value\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_median,[-180,180],[-90,90] );

end
%% mean year (satellite)
clc
clear
for year=2015:2023


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']);

    NEE_mean(:,:,1)=A_NEE;
    NEE_mean(:,:,2)=B_NEE;
    NEE_mean(:,:,3)=C_NEE;
    NEE_mean(:,:,4)=D_NEE;

    NEE_mean=nanmean(NEE_mean,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_mean,[-180,180],[-90,90] );


end
%% mean year (satellite_plus)
clc
clear
for year=2015:2023


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\year\CAMSS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\year\CMSF_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\year\GCASv2_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\year\GONGGA_NEE_' num2str(year) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\year\COLA_NEE_' num2str(year) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\year\NTFVAR_NEE_' num2str(year) '.tif']);

    NEE_mean(:,:,1)=A_NEE;
    NEE_mean(:,:,2)=B_NEE;
    NEE_mean(:,:,3)=C_NEE;
    NEE_mean(:,:,4)=D_NEE;
    NEE_mean(:,:,5)=E_NEE;
    NEE_mean(:,:,6)=F_NEE;

    NEE_mean=nanmean(NEE_mean,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_mean,[-180,180],[-90,90] );


end
%% mean year (non-satellite)

for year=2015:2023


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\year\CAMS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\year\CarboScope_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\year\COLA_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\year\CTE_NEE_' num2str(year) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\year\CT_NEE_' num2str(year) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\year\IAPCAS_NEE_' num2str(year) '.tif']);
    G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\year\MIROC-ACTM_NEE_' num2str(year) '.tif']);
    H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\year\NISMON_NEE_' num2str(year) '.tif']);
    I_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\year\NTFVAR_NEE_' num2str(year) '.tif']);
    J_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\year\UoE_NEE_' num2str(year) '.tif']);



    NEE_mean(:,:,1)=A_NEE;
    NEE_mean(:,:,2)=B_NEE;
    NEE_mean(:,:,3)=C_NEE;
    NEE_mean(:,:,4)=D_NEE;
    NEE_mean(:,:,5)=E_NEE;
    NEE_mean(:,:,6)=F_NEE;
    NEE_mean(:,:,7)=G_NEE;
    NEE_mean(:,:,8)=H_NEE;
    NEE_mean(:,:,9)=I_NEE;
    NEE_mean(:,:,10)=J_NEE;

    NEE_mean=nanmean(NEE_mean,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_mean,[-180,180],[-90,90] );

end
%% mean year (non-satellite_PLUS)
clc
clear
for year=2015:2023


    A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\year\CAMS_NEE_' num2str(year) '.tif']);
    B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\year\CarboScope_NEE_' num2str(year) '.tif']);
    C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\year\CTE_NEE_' num2str(year) '.tif']);
    D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\year\CT_NEE_' num2str(year) '.tif']);
    E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\year\IAPCAS_NEE_' num2str(year) '.tif']);
    F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\year\MIROC-ACTM_NEE_' num2str(year) '.tif']);
    G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\year\NISMON_NEE_' num2str(year) '.tif']);
    H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\year\UoE_NEE_' num2str(year) '.tif']);

    NEE_mean(:,:,1)=A_NEE;
    NEE_mean(:,:,2)=B_NEE;
    NEE_mean(:,:,3)=C_NEE;
    NEE_mean(:,:,4)=D_NEE;
    NEE_mean(:,:,5)=E_NEE;
    NEE_mean(:,:,6)=F_NEE;
    NEE_mean(:,:,7)=G_NEE;
    NEE_mean(:,:,8)=H_NEE;

    NEE_mean=nanmean(NEE_mean,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\year\NEE_' num2str(year) '.tif']
    wtif(filename,NEE_mean,[-180,180],[-90,90] );

end
%% mean mmonth
clc
clear

for year=2015:2023

    for month=1:12
        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year) '_' num2str(month) '.tif']);
        G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year) '_' num2str(month) '.tif']);
        H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        I_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        J_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        K_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year) '_' num2str(month) '.tif']);
        L_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year) '_' num2str(month) '.tif']);
        M_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year) '_' num2str(month) '.tif']);
        N_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year) '_' num2str(month) '.tif']);



        NEE_mean(:,:,1)=A_NEE;
        NEE_mean(:,:,2)=B_NEE;
        NEE_mean(:,:,3)=C_NEE;
        NEE_mean(:,:,4)=D_NEE;
        NEE_mean(:,:,5)=E_NEE;
        NEE_mean(:,:,6)=F_NEE;
        NEE_mean(:,:,7)=G_NEE;
        NEE_mean(:,:,8)=H_NEE;
        NEE_mean(:,:,9)=I_NEE;
        NEE_mean(:,:,10)=J_NEE;
        NEE_mean(:,:,11)=K_NEE;
        NEE_mean(:,:,12)=L_NEE;
        NEE_mean(:,:,13)=M_NEE;
        NEE_mean(:,:,14)=N_NEE;


        NEE_median=nanmedian(NEE_mean,3);
        NEE_mean=nanmean(NEE_mean,3);

        filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,NEE_mean,[-180,180],[-90,90] );
        filename=['E:\phd_file\Boreal_North_America\GCB2024\Median_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,NEE_median,[-180,180],[-90,90] );


    end
end
%% mean month (satellite)
clc
clear

for year=2015:2023

    for month=1:12
        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);




        NEE_mean(:,:,1)=A_NEE;
        NEE_mean(:,:,2)=B_NEE;
        NEE_mean(:,:,3)=C_NEE;
        NEE_mean(:,:,4)=D_NEE;

        NEE_mean=nanmean(NEE_mean,3);

        filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,NEE_mean,[-180,180],[-90,90] );


    end
end
%% mean month (satellite_PLUS)
clc
clear

for year=2015:2023

    for month=1:12
        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS-satellite\month\CAMSS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CMS-Flux\month\CMSF_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GCASv2\month\GCASv2_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\GONGGA\month\GONGGA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\COLA\month\COLA_NEE_' num2str(year) '_' num2str(month) '.tif']);
        F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NTFVAR\month\NTFVAR_NEE_' num2str(year) '_' num2str(month) '.tif']);


        NEE_mean(:,:,1)=A_NEE;
        NEE_mean(:,:,2)=B_NEE;
        NEE_mean(:,:,3)=C_NEE;
        NEE_mean(:,:,4)=D_NEE;
        NEE_mean(:,:,5)=E_NEE;
        NEE_mean(:,:,6)=F_NEE;


        NEE_mean=nanmean(NEE_mean,3);

        filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,NEE_mean,[-180,180],[-90,90] );


    end
end
%% mean month (non-satellite-PLUS)
clc
clear

for year=2015:2023

    for month=1:12
        A_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CAMS\month\CAMS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        B_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CarbonScope\month\CarboScope_NEE_' num2str(year) '_' num2str(month) '.tif']);
        C_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CTE\month\CTE_NEE_' num2str(year) '_' num2str(month) '.tif']);
        D_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\CT-NOAA\month\CT_NEE_' num2str(year) '_' num2str(month) '.tif']);
        E_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\IAPCAS\month\IAPCAS_NEE_' num2str(year) '_' num2str(month) '.tif']);
        F_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\MIROC-ACTM\month\MIROC-ACTM_NEE_' num2str(year) '_' num2str(month) '.tif']);
        G_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\NISMON-CO2\month\NISMON_NEE_' num2str(year) '_' num2str(month) '.tif']);
        H_NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\UoE\month\UoE_NEE_' num2str(year) '_' num2str(month) '.tif']);



        NEE_mean(:,:,1)=A_NEE;
        NEE_mean(:,:,2)=B_NEE;
        NEE_mean(:,:,3)=C_NEE;
        NEE_mean(:,:,4)=D_NEE;
        NEE_mean(:,:,5)=E_NEE;
        NEE_mean(:,:,6)=F_NEE;
        NEE_mean(:,:,7)=G_NEE;
        NEE_mean(:,:,8)=H_NEE;


        NEE_mean=nanmean(NEE_mean,3);

        filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,NEE_mean,[-180,180],[-90,90] );

    end
end
%% mean summer satellite (6-8)
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\summer\NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% mean summer satellite-PLUS (6-8)
clc
clear

for year=2015:2023

    for month=6:8
        NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        Data_list(:,:,month-5)=NEE_temp;
    end
    result=nansum(Data_list,3);

    filename=['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\summer\NEE_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear Data_list result
end
%% mean ER year
clc
clear

for year=2015:2023

    NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\year\NEE_' num2str(year) '.tif']);
    GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    result=NEE+GPP;
    % result(result<0)=nan;
    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\mean\year\ER_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear  result
end
%% mean ER year (satellite)
clc
clear

for year=2015:2023

    NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\year\NEE_' num2str(year) '.tif']);
    GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    result=NEE+GPP;
    % result(result<0)=nan;
    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\year\ER_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear  result
end
%% mean ER year (satellite-plus)
clc
clear

for year=2015:2023

    NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\year\NEE_' num2str(year) '.tif']);
    GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    result=NEE+GPP;
    % result(result<0)=nan;
    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\year\ER_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear  result
end
%% mean ER year (non-satellite)
clc
clear

for year=2015:2023

    NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite\year\NEE_' num2str(year) '.tif']);
    GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    result=NEE+GPP;
    % result(result<0)=nan;
    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite\year\ER_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear  result
end
%% mean ER year (non-satellite-plus)
clc
clear

for year=2015:2023

    NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\year\NEE_' num2str(year) '.tif']);
    GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\year\GPP_' num2str(year) '.tif']);

    result=NEE+GPP;
    % result(result<0)=nan;
    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\year\ER_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear  result
end
%% mean ER mmonth
clc
clear

for year=2015:2023

    for month=1:12
        NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);

        result=NEE+GPP;
        % result(result<0)=nan;

        filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\mean\month\ER_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,result,[-180,180],[-90,90] );
        clear  result
    end
end
%% mean ER month (satellite)
clc
clear

for year=2015:2023

    for month=1:12
        NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);

        result=NEE+GPP;
        % result(result<0)=nan;

        filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\month\ER_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,result,[-180,180],[-90,90] );
        clear  result
    end
end
%% mean ER month (satellite-plus)
clc
clear

for year=2015:2023

    for month=1:12
        NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);

        result=NEE+GPP;
        % result(result<0)=nan;

        filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,result,[-180,180],[-90,90] );
        clear  result
    end
end
%% mean ER month (non-satellite)
clc
clear

for year=2015:2023

    for month=1:12
        NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);

        result=NEE+GPP;
        % result(result<0)=nan;

        filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite\month\ER_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,result,[-180,180],[-90,90] );
        clear  result
    end
end
%% mean ER month (non-satellite-plus)
clc
clear

for year=2015:2023

    for month=1:12
        NEE=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\non_satellite_plus\month\NEE_' num2str(year) '_' num2str(month) '.tif']);
        GPP=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\month\GPP_' num2str(year) '_' num2str(month) '.tif']);

        result=NEE+GPP;
        % result(result<0)=nan;

        filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\non-satellite_plus\month\ER_' num2str(year) '_' num2str(month) '.tif']
        wtif(filename,result,[-180,180],[-90,90] );
        clear  result
    end
end
%% mean ER summer satellite (6-8)
clc
clear

for year=2015:2023


    NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite\summer\NEE_summer_' num2str(year) '.tif']);
    GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_' num2str(year) '.tif']);

    result=NEE_temp+GPP_temp;

    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite\summer\ER_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear result
end
%% mean ER summer satellite-plus (6-8)
clc
clear

for year=2015:2023


    NEE_temp=importdata(['E:\phd_file\Boreal_North_America\GCB2024\Mean_value\satellite_plus\summer\NEE_summer_' num2str(year) '.tif']);
    GPP_temp=importdata(['E:\phd_file\Boreal_North_America\GPP\mean_value\summer\GPP_summer_' num2str(year) '.tif']);

    result=NEE_temp+GPP_temp;

    filename=['E:\phd_file\Boreal_North_America\GCB2024\ER\satellite_plus\summer\ER_summer_' num2str(year) '.tif']
    wtif(filename,result,[-180,180],[-90,90] );
    clear result
end