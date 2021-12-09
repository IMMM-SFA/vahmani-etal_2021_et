clear 
close all
clc
% addpath('/Volumes/GoogleDrive/My Drive/5.Projects/7.WRF-UCM_WaterDemand_ClimateChange/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/5.Projects/6.LBNL_Bay_Area_Climate_Readiness_Project_part3/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/5.Projects/2016_Ultra_High_Res_Climate_Modeling_LDRD_and_MLA_part3/2.WRF_PostProcessing/0.Pouya_matlab_fx')
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

%loop over variables 
Variables={'Irr'   ,'ET'        ,'RnStar'    ,'D'  ,'s'    ,'ra' ,'Ta','rs' };
Units    ={'mm/day','mm/12hrDay','mm/12hrDay','kPa','kPa/C','s/m','C' ,'s/m'};

for var_i=1:length(Variables)
varname=Variables{var_i};
varunit=Units    {var_i};
eval(sprintf('base=%s1;est=%s2;',varname,varname));

% statistics
stat_base= pouya_time_series2( base, data_LU_INDEX, LULC, LULCl );
stat_est = pouya_time_series2( est , data_LU_INDEX, LULC, LULCl );

%calc changes
delta=(est-base);
% statistics
stat_delta= pouya_time_series2( delta, data_LU_INDEX, LULC, LULCl );

% plot spatial options
urban_border=1;
if strcmp(varname,'Irr')
urban_border=0;
end
county_border=0;
state_border=1;

% plot spatial base
plot_data=nanmean(base,3);
plot_data(data_LU_INDEX==17)=NaN;
cmin_diff=prctile(plot_data(:),1 );
cmax_diff=prctile(plot_data(:),99);
cpt=NaN;%'YlOrRd_09';
pouya_fig_geoshow_with_cptcmapOption_agBorder(cpt,plot_data,data_XLAT,data_XLONG,data_LU_INDEX,cmin_diff,cmax_diff,urban_border,county_border,state_border);
title(sprintf('%s (%s) \n%f(%s),%f(%s),%f(%s),%f(%s)',varname,varunit,nanmean(stat_base.agric),LULCl{1,1},nanmean(stat_base.urb),LULCl{2,1},nanmean(stat_base.land_nonurb_nonag),LULCl{3,1},nanmean(stat_base.water),LULCl{4,1}),'FontSize',9);
savename=sprintf('%s/%s_Base_Spatial',savefolder,varname);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

% plot spatial changes red blue
plot_data=nanmean(delta,3);
plot_data(data_LU_INDEX==17)=NaN;
plot_data(plot_data==0)=NaN;%make zeros white
cmin_diff=-1*max(abs(prctile(plot_data(:),0.1)),abs(prctile(plot_data(:),99.9)));
cmax_diff=   max(abs(prctile(plot_data(:),0.1)),abs(prctile(plot_data(:),99.9)));
cpt='Oranges_09_strong_Blues_09';
pouya_fig_geoshow_with_cptcmapOption_agBorder(cpt,plot_data,data_XLAT,data_XLONG,data_LU_INDEX,cmin_diff,cmax_diff,urban_border,county_border,state_border);
title(sprintf('%s (%s) change\n%f(%s),%f(%s),%f(%s),%f(%s)',varname,varunit,nanmean(stat_delta.agric),LULCl{1,1},nanmean(stat_delta.urb),LULCl{2,1},nanmean(stat_delta.land_nonurb_nonag),LULCl{3,1},nanmean(stat_delta.water),LULCl{4,1}),'FontSize',9);
savename=sprintf('%s/%s_Change_Spatial',savefolder,varname);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

%calc percent changes
delta_perc=(est-base)./base*100;
%avoid unrealistic percentage changes over close to zero baselines (this mostly affect error bars)
delta_perc(base<=prctile(base(:),0.1))=NaN;
fprintf('note: %s < %f (%s) were not condidered in perc change calc\n',varname,prctile(base(:),0.1),varunit);
%for rs added another layer cause a few values very super high in rs2
if strcmp(varname,'rs')
delta_perc(est>prctile(est(:),99.9))=NaN;
fprintf('note: %s > %f (%s) were not condidered in perc change calc\n',varname,prctile(est(:),99.9),varunit);
end

% statistics
stat_delta_perc= pouya_time_series2( delta_perc,data_LU_INDEX, LULC, LULCl );

% plot spatial percent change
plot_data=nanmean(delta_perc,3);
plot_data(data_LU_INDEX==17)=NaN;
plot_data(plot_data==0)=NaN;%make zeros white
cmin_diff=-30;
cmax_diff= 30;
cpt='Oranges_09_strong_Blues_09';
pouya_fig_geoshow_with_cptcmapOption_agBorder(cpt,plot_data,data_XLAT,data_XLONG,data_LU_INDEX,cmin_diff,cmax_diff,urban_border,county_border,state_border);
title(sprintf('%s percent change\n%2.0f(%s),%2.0f(%s),%2.0f(%s),%2.0f(%s)',varname,nanmean(stat_delta_perc.agric),LULCl{1,1},nanmean(stat_delta_perc.urb),LULCl{2,1},nanmean(stat_delta_perc.land_nonurb_nonag),LULCl{3,1},nanmean(stat_delta_perc.water),LULCl{4,1}),'FontSize',9);
savename=sprintf('%s/%s_PercChange_Spatial',savefolder,varname);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

% plot bar percent change
plot_data    =[nanmean(stat_delta_perc.agric) nanmean(stat_delta_perc.urb) nanmean(stat_delta_perc.land_nonurb_nonag)]; 
plot_data_err=[nanstd( stat_delta_perc.agric) nanstd( stat_delta_perc.urb) nanstd( stat_delta_perc.land_nonurb_nonag)]; 
figure('position',[1000,1000,250,500]);
bar(1:3,plot_data,'FaceColor',[0.9290,0.6940,0.1250])   
hold on
er = errorbar(1:3,plot_data,plot_data_err,plot_data_err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ylim([-5 45]);
set(gca,'YTick',(-10:10:40));
set(gca,'YTickLabel',{'-10%','0','10%','20%','30%','40%'});
set(gca,'YAxisLocation','right');
set(gca,'XTickLabel',{'ag.','urb.','land'});
set(gca,'FontSize',18);
hold off
title(sprintf('%s percent change\n%2.0f(%s),%2.0f(%s),%2.0f(%s),%2.0f(%s)',varname,nanmean(stat_delta_perc.agric),LULCl{1,1},nanmean(stat_delta_perc.urb),LULCl{2,1},nanmean(stat_delta_perc.land_nonurb_nonag),LULCl{3,1},nanmean(stat_delta_perc.water),LULCl{4,1}),'FontSize',9);
savename=sprintf('%s/%s_PercChange_Bar',savefolder,varname);
saveas(gca,strcat(savename,'.fig'));
set(gcf,'Units','pixels');scrpos = get(gcf,'Position');newpos = scrpos/100;set(gcf,'PaperUnits','inches','PaperPosition',newpos);%Make the output file the same size as the figure on the screen
print('-dpng', strcat(savename,'.png'), '-r300');

%save data
savename=sprintf('%s/%s_StatData.mat',savefolder,varname);
save(savename,'varname','varunit','stat_base','stat_est','stat_delta','stat_delta_perc','-v7.3')
end

end