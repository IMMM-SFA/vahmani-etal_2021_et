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

% plot bar percent change
plot_data    =[nanmean(stat_tot_ET    .agric) nanmean(stat_tot_ET    .urb) nanmean(stat_tot_ET    .land_nonurb_nonag);...
               nanmean(stat_att_RnStar.agric) nanmean(stat_att_RnStar.urb) nanmean(stat_att_RnStar.land_nonurb_nonag);...
               nanmean(stat_att_D     .agric) nanmean(stat_att_D     .urb) nanmean(stat_att_D     .land_nonurb_nonag);...
               nanmean(stat_att_rs    .agric) nanmean(stat_att_rs    .urb) nanmean(stat_att_rs    .land_nonurb_nonag);...
               nanmean(stat_att_ra    .agric) nanmean(stat_att_ra    .urb) nanmean(stat_att_ra    .land_nonurb_nonag);...
               nanmean(stat_att_s     .agric) nanmean(stat_att_s     .urb) nanmean(stat_att_s     .land_nonurb_nonag)]; 
plot_data_err=[nanstd( stat_tot_ET    .agric) nanstd( stat_tot_ET    .urb) nanstd( stat_tot_ET    .land_nonurb_nonag);...
               nanstd( stat_att_RnStar.agric) nanstd( stat_att_RnStar.urb) nanstd( stat_att_RnStar.land_nonurb_nonag);...
               nanstd( stat_att_D     .agric) nanstd( stat_att_D     .urb) nanstd( stat_att_D     .land_nonurb_nonag);...
               nanstd( stat_att_rs    .agric) nanstd( stat_att_rs    .urb) nanstd( stat_att_rs    .land_nonurb_nonag);...
               nanstd( stat_att_ra    .agric) nanstd( stat_att_ra    .urb) nanstd( stat_att_ra    .land_nonurb_nonag);...
               nanstd( stat_att_s     .agric) nanstd( stat_att_s     .urb) nanstd( stat_att_s     .land_nonurb_nonag)];
           
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
ylim([-0.4 0.4]);
set(gca,'XTickLabel',{'ET','RnStar_a_t_t','D_a_t_t','rs_a_t_t','ra_a_t_t','s_a_t_t'});
set(gca,'FontSize',12);
title(sprintf('ET change and att (mm/day)\n%s',filename{1}),'FontSize',10);
savename=sprintf('%s/All_Att_Bar',savefolder);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

%save data
savename=sprintf('%s/All_Att_StatData.mat',savefolder);
save(savename,'savename','plot_data','plot_data_err','-v7.3')

end