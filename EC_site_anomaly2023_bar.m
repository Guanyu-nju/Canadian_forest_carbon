clc
clear

% SOB
SOB_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\SOB_monthly.xlsx"); clear year
SOB_data.Year = year(SOB_data.Time_UTC);
SOB_annual_sums = varfun(@sum, SOB_data, 'GroupingVariables', 'Year', ...
                     'InputVariables', {'EcosystemRespiration','GrossEcosystemPhotosynthesis', 'NetEcosystemProduction'});
SOB_NEP_2023anomaly=SOB_annual_sums{end,5}-nanmean(SOB_annual_sums{1:end-1,5});
SOB_GPP_2023anomaly=SOB_annual_sums{end,4}-nanmean(SOB_annual_sums{1:end-1,4});
SOB_ER_2023anomaly=SOB_annual_sums{end,3}-nanmean(SOB_annual_sums{1:end-1,3});
SOB_result=[SOB_NEP_2023anomaly,SOB_GPP_2023anomaly,SOB_ER_2023anomaly];

% OJP
OJP_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\OJP_monthly.xlsx"); clear year
OJP_data.Year = year(OJP_data.Time_UTC);
OJP_annual_sums = varfun(@sum, OJP_data, 'GroupingVariables', 'Year', ...
                     'InputVariables', {'EcosystemRespiration','GrossEcosystemPhotosynthesis', 'NetEcosystemProduction'});
OJP_NEP_2023anomaly=OJP_annual_sums{end,5}-nanmean(OJP_annual_sums{1:end-1,5});
OJP_GPP_2023anomaly=OJP_annual_sums{end,4}-nanmean(OJP_annual_sums{1:end-1,4});
OJP_ER_2023anomaly=OJP_annual_sums{end,3}-nanmean(OJP_annual_sums{1:end-1,3});
OJP_result=[OJP_NEP_2023anomaly,OJP_GPP_2023anomaly,OJP_ER_2023anomaly];

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
SCC_result=[SCC_NEP_2023anomaly,SCC_GPP_2023anomaly,SCC_ER_2023anomaly];

% HPC
HPC_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\Guanyu_ForestSites_2015-2023.xlsx",'Sheet', 'HPC');
HPC_annual_sums = varfun(@sum, HPC_data, 'GroupingVariables', 'year', ...
                     'InputVariables', {'nee','gpp', 'reco'});

HPC_annual_sums{:,3}=-HPC_annual_sums{:,3};
HPC_annual_sums{:,4}=HPC_annual_sums{:,3}+HPC_annual_sums{:,5};
HPC_NEP_2023anomaly=HPC_annual_sums{end,3}-nanmean(HPC_annual_sums{1:end-1,3});
HPC_GPP_2023anomaly=HPC_annual_sums{end,4}-nanmean(HPC_annual_sums{1:end-1,4});
HPC_ER_2023anomaly=HPC_annual_sums{end,5}-nanmean(HPC_annual_sums{1:end-1,5});
HPC_result=[HPC_NEP_2023anomaly,HPC_GPP_2023anomaly,HPC_ER_2023anomaly];

% TPD
for year=2015:2023
    data_temp=readtable(['E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\CA-TPD\TPD_Monthly_' num2str(year) '.csv']);
    TPD_annual_sum(year-2014,1)=sum(data_temp{1:12,4})+sum(data_temp{1:12,5});
    TPD_annual_sum(year-2014,2)=sum(data_temp{1:12,4});
    TPD_annual_sum(year-2014,3)=sum(data_temp{1:12,5});

end
TPD_NEP_2023anomaly=TPD_annual_sum(end,3)-mean(TPD_annual_sum(1:8,3))
TPD_GPP_2023anomaly=TPD_annual_sum(end,1)-mean(TPD_annual_sum(1:8,1))
TPD_ER_2023anomaly=TPD_annual_sum(end,2)-mean(TPD_annual_sum(1:8,2))
TPD_result=[TPD_NEP_2023anomaly,TPD_GPP_2023anomaly,TPD_ER_2023anomaly];

% TP3
for year=2015:2023
    data_temp=readtable(['E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\CA-TPD\CA_TP3_Monthly_' num2str(year) '.csv']);
    TP3_annual_sum(year-2014,1)=sum(data_temp{1:12,4})+sum(data_temp{1:12,5});
    TP3_annual_sum(year-2014,2)=sum(data_temp{1:12,4});
    TP3_annual_sum(year-2014,3)=sum(data_temp{1:12,5});

end
TP3_NEP_2023anomaly=TP3_annual_sum(end,3)-mean(TP3_annual_sum(1:8,3))
TP3_GPP_2023anomaly=TP3_annual_sum(end,1)-mean(TP3_annual_sum(1:8,1))
TP3_ER_2023anomaly=TP3_annual_sum(end,2)-mean(TP3_annual_sum(1:8,2))
TP3_result=[TP3_NEP_2023anomaly,TP3_GPP_2023anomaly,TP3_ER_2023anomaly];

% 
% 
SOB_data=readtable("E:\phd_file\Boreal_North_America\Ground_flux\EC_Flux\Cbo_monthly.xlsx"); clear year
SOB_data.Year = year(SOB_data.Date_Time);
SOB_annual_sums = varfun(@sum, SOB_data, 'GroupingVariables', 'Year', ...
                     'InputVariables', {'NEE','ER', 'GPP'});
SOB_NEP_2023anomaly=SOB_annual_sums{end,3}-nanmean(SOB_annual_sums{1:end-1,3});
SOB_GPP_2023anomaly=SOB_annual_sums{end,5}-nanmean(SOB_annual_sums{1:end-1,5});
SOB_ER_2023anomaly=SOB_annual_sums{end,4}-nanmean(SOB_annual_sums{1:end-1,4});
SOB_result=[SOB_NEP_2023anomaly,SOB_GPP_2023anomaly,SOB_ER_2023anomaly];

EC_result=[HPC_result;SCC_result;SOB_result;TP3_result;TPD_result];

%%
f=figure
set(gcf,'unit','pixels','position',[836,866,630,360]);
% t = tiledlayout(1,1);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
b=bar(EC_result,'FaceColor','flat'); hold on
b(1).CData = [44,107,179]/255;
b(2).CData = [17,119,51]/255;
b(3).CData = [193,0,1]/255;
ylim([-150,150])
set(gca,'YTick', [-100:50:150]);
ylabel('Anomaly (gC m^{-2} yr^{-1})','FontName','Arial','FontSize',14);
% xlabel('Site name','FontName','Arial','FontSize',14)
% set(gca,'XTick',[1:1:3],'FontName','Arial','fontsize',14)
set(gca,'XTickLabel',{'CA-HPC','CA-SCC','CA-SOB','CA-TP3','CA-TPD'},'FontName','Arial','fontsize',14)
legend({'NEP','GPP','TER'},'NumColumns',3,'FontName','Arial','FontSize',14,'Box','on','Location','northeast')
% text('string','a','Units','normalized','position',[-0.158500470010452 1.01781870738078 0],'FontName','Arial','FontSize',20,'fontweight','bold')

result=['E:\phd_file\Boreal_North_America\Result\V9\EC_site_anomaly2023_bar.png']
% print(result,f,'-r600','-dpng');
