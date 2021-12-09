clear 
close all
clc
% addpath('/Volumes/GoogleDrive/My Drive/5.Projects/7.WRF-UCM_WaterDemand_ClimateChange/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/5.Projects/6.LBNL_Bay_Area_Climate_Readiness_Project_part3/2.WRF_PostProcessing/0.Pouya_matlab_fx')
%%
for filename={'Data_partial_differentials_CNRMcoolRoof-CNRM_v11',...
              'Data_partial_differentials_HADGcoolRoof-HADG_v11',...
              'Data_partial_differentials_HADG-baseline_v10',...
              'Data_partial_differentials_CNRM-baseline_v10',...
              'Data_partial_differentials_HADGcoolRoof-baseline_v10',...
              'Data_partial_differentials_CNRMcoolRoof-baseline_v10'};
      
load(filename{1});
savefolder=sprintf('%s/%s',savefolder,filename{1});
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end
%%
%calc perc change and avoid unrealistic percentage changes over close to zero baselines (this mostly affect error bars)
deltaET_perc    =(ET2    -ET1    )./ET1    *100; deltaET_perc    (ET1    <prctile(ET1    (:),0.1))=NaN;
deltaRnStar_perc=(RnStar2-RnStar1)./RnStar1*100; deltaRnStar_perc(RnStar1<prctile(RnStar1(:),0.1))=NaN;
deltaD_perc     =(D2     -D1     )./D1     *100; deltaD_perc     (D1     <prctile(D1     (:),0.1))=NaN;
deltara_perc    =(ra2    -ra1    )./ra1    *100; deltara_perc    (ra1    <prctile(ra1    (:),0.1))=NaN;
deltas_perc     =(s2     -s1     )./s1     *100; deltas_perc     (s1     <prctile(s1     (:),0.1))=NaN;
deltars_perc    =(rs2    -rs1    )./rs1    *100; deltars_perc    (rs1    <prctile(rs1    (:),0.1))=NaN;
%for rs added another layer cause a few values very super high in rs2
deltars_perc    (rs2    >prctile(rs2    (:),99.9))=NaN;

% statistics
stat_change_ET    = pouya_time_series2( deltaET_perc    ,data_LU_INDEX, LULC, LULCl );
stat_change_RnStar= pouya_time_series2( deltaRnStar_perc,data_LU_INDEX, LULC, LULCl );
stat_change_D     = pouya_time_series2( deltaD_perc     ,data_LU_INDEX, LULC, LULCl );
stat_change_rs    = pouya_time_series2( deltars_perc    ,data_LU_INDEX, LULC, LULCl );
stat_change_ra    = pouya_time_series2( deltara_perc    ,data_LU_INDEX, LULC, LULCl );
stat_change_s     = pouya_time_series2( deltas_perc     ,data_LU_INDEX, LULC, LULCl );

% plot bar percent change
plot_data    =[nanmean(stat_change_ET    .agric) nanmean(stat_change_ET    .urb) nanmean(stat_change_ET    .land_nonurb_nonag);...
               nanmean(stat_change_RnStar.agric) nanmean(stat_change_RnStar.urb) nanmean(stat_change_RnStar.land_nonurb_nonag);...
               nanmean(stat_change_D     .agric) nanmean(stat_change_D     .urb) nanmean(stat_change_D     .land_nonurb_nonag);...
               nanmean(stat_change_rs    .agric) nanmean(stat_change_rs    .urb) nanmean(stat_change_rs    .land_nonurb_nonag);...
               nanmean(stat_change_ra    .agric) nanmean(stat_change_ra    .urb) nanmean(stat_change_ra    .land_nonurb_nonag);...
               nanmean(stat_change_s     .agric) nanmean(stat_change_s     .urb) nanmean(stat_change_s     .land_nonurb_nonag)]; 

plot_data_err=[nanstd( stat_change_ET    .agric) nanstd( stat_change_ET    .urb) nanstd( stat_change_ET    .land_nonurb_nonag);...
               nanstd( stat_change_RnStar.agric) nanstd( stat_change_RnStar.urb) nanstd( stat_change_RnStar.land_nonurb_nonag);...
               nanstd( stat_change_D     .agric) nanstd( stat_change_D     .urb) nanstd( stat_change_D     .land_nonurb_nonag);...
               nanstd( stat_change_rs    .agric) nanstd( stat_change_rs    .urb) nanstd( stat_change_rs    .land_nonurb_nonag);...
               nanstd( stat_change_ra    .agric) nanstd( stat_change_ra    .urb) nanstd( stat_change_ra    .land_nonurb_nonag);...
               nanstd( stat_change_s     .agric) nanstd( stat_change_s     .urb) nanstd( stat_change_s     .land_nonurb_nonag)];
           
figure;
b=bar(plot_data);
b(1).FaceColor=[0.8500,0.3250,0.0980];
b(2).FaceColor=[0.0000,0.4470,0.7410];
b(3).FaceColor=[0.9290,0.6940,0.1250];
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(plot_data, 1);
nbars = size(plot_data, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, plot_data(:,i), plot_data_err(:,i), 'k', 'linestyle', 'none');
end
hold off
ylim([-30 30]);
set(gca,'YTick',(-30:10:30));
set(gca,'YTickLabel',{'-30%','-20%','-10%','0','10%','20%','30%','40%'});
set(gca,'XTickLabel',{'ET','RnStar_a_t_t','D_a_t_t','rs_a_t_t','ra_a_t_t','s_a_t_t'});
set(gca,'FontSize',12);
title(sprintf('ET change and att (mm/day)\n%s',filename{1}),'FontSize',10);
savename=sprintf('%s/All_PercChange_Bar',savefolder);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

%save data
savename=sprintf('%s/All_PercChange_StatData.mat',savefolder);
save(savename,'savename','plot_data','plot_data_err','-v7.3')
end