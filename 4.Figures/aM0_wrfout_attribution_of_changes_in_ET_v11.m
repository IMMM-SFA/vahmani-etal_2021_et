clear 
close all
clc
% addpath('/Volumes/GoogleDrive/My Drive/5.Projects/7.WRF-UCM_WaterDemand_ClimateChange/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/5.Projects/6.LBNL_Bay_Area_Climate_Readiness_Project_part3/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/5.Projects/2016_Ultra_High_Res_Climate_Modeling_LDRD_and_MLA_part3/2.WRF_PostProcessing/0.Pouya_matlab_fx')
addpath('/Volumes/GoogleDrive/My Drive/PROJECTS/2016_Ultra_High_Res_Climate_Modeling_LDRD_and_MLA/STUDY 2018 Climate Change Extreme Heat Exposure CA/2.WRF_PostProcessing/0.Pouya_matlab_fx')
%%
% folder that processed data is saved
runname1='S31_base_ens01_WRFRun_';plot_runname1='BASE'        ;pop_scenario1='Y2010A2';
runname2='S31_cnrm_ens01_WRFRun_';plot_runname2='CNRM'        ;pop_scenario2='Y2010A2';
runname3='S41_cnrm_ens01_WRFRun_';plot_runname3='CNRMcoolRoof';pop_scenario3='Y2010A2';
runname4='S31_hadg_ens01_WRFRun_';plot_runname4='HADG'        ;pop_scenario4='Y2010A2';
runname5='S41_hadg_ens01_WRFRun_';plot_runname5='HADGcoolRoof';pop_scenario5='Y2010A2';

foldername='/Volumes/PVahmani_G_RAID/Projects/6.LBNL_Bay_Area_Climate_Readiness_Project_part3/2.WRF_Outputs_daily_daytime';
savefolder='/Volumes/GoogleDrive/My Drive/5.Projects/7.WRF-UCM_WaterDemand_ClimateChange/2.WRF_PostProcessing/2.Figures/3.ETattribution';
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%land cover land use of interest
hiUrb=24;meUrb=25;loUrb=26;
agric=38;
water=17;
land=100;land_nonurb=101;land_nonag=102;land_nonurb_nonag=103;
all=1000;
%set the land cover land use of interest for each domain
LULC ={agric;[hiUrb;meUrb;loUrb];land_nonurb_nonag;water};clear hiUrb meUrb loUrb agric water land land_nonurb land_nonag land_nonurb_nonag all
LULCc={'g';'r';'k';'b'};
LULCl={'agric';'urb';'land_nonurb_nonag';'water'};
%% load basics
%data_XLAT,data_XLONG,data_LU_INDEX
[ ~,~,~,~,data_XLAT,data_XLONG,data_LU_INDEX ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname1,'ET2'           );
DATA_LU_INDEX=zeros(size(data_LU_INDEX,1),size(data_LU_INDEX,2),15);
for k=1:15
    DATA_LU_INDEX(:,:,k)=data_LU_INDEX;
end;clear k
DATA_LU_INDEX=single(DATA_LU_INDEX);
%Urban Fraction (-)
[ data_FRC_URB2D,~,~,~,~,~,~] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname1,'FRC_URB2D'     );
data_FRC_URB2D=data_FRC_URB2D(:,:,1,1);
DATA_FRC_URB2D=zeros(size(data_LU_INDEX,1),size(data_LU_INDEX,2),15);
for k=1:15
    DATA_FRC_URB2D(:,:,k)=data_FRC_URB2D;
end;clear k
DATA_FRC_URB2D=single(DATA_FRC_URB2D);
%only CA and ocean 
ca_poly_on_in_withOcean = pouya_CA_only_tight_withOcean( data_XLONG,data_XLAT );
CA_poly_on_in_withOcean=zeros(size(ca_poly_on_in_withOcean,1),size(ca_poly_on_in_withOcean,2),15);
for k=1:15
    CA_poly_on_in_withOcean(:,:,k)=ca_poly_on_in_withOcean;
end;clear k ca_poly_on_in_withOcean
CA_poly_on_in_withOcean=single(CA_poly_on_in_withOcean);
%% baseline
runname_1=runname2;
plot_runname_1=plot_runname2;
%TSK (C) (calc and use T1 for urban which reflect veg (pervious) portion of urban gridcells
[ data1_in_daily_TSK          ,~,varunit_TSK   ,vardesc_TSK   ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'TSK'          );
data1_in_daily_TSK =squeeze(nanmean(data1_in_daily_TSK ,3));
[ data1_in_daily_TS_URB       ,~,~             ,~             ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'TS_URB'       );
data1_in_daily_TS_URB =squeeze(nanmean(data1_in_daily_TS_URB ,3));
data1_in_daily_T1=(data1_in_daily_TSK-data1_in_daily_TS_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data1_in_daily_TS_URB
data1_in_daily_TSK(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data1_in_daily_T1(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data1_in_daily_T1

%HFX (W/m2 or J/m2.s) (calc and use SHEAT for urban which reflect veg (pervious) portion of urban gridcells
[ data1_in_daily_HFX          ,~,varunit_HFX   ,vardesc_HFX   ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'HFX'          );
data1_in_daily_HFX =squeeze(nanmean(data1_in_daily_HFX ,3));
[ data1_in_daily_SH_URB       ,~,varunit_SH_URB,vardesc_SH_URB,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'SH_URB'       );
data1_in_daily_SH_URB =squeeze(nanmean(data1_in_daily_SH_URB ,3));
data1_in_daily_SHEAT=(data1_in_daily_HFX-data1_in_daily_SH_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data1_in_daily_SH_URB
data1_in_daily_HFX(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data1_in_daily_SHEAT(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data1_in_daily_SHEAT

%T2 (C)
[ data1_in_daily_T2           ,~,varunit_T2    ,vardesc_T2    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'T2'           );
data1_in_daily_T2 =squeeze(nanmean(data1_in_daily_T2 ,3));
data1_in_daily_T2           (~CA_poly_on_in_withOcean)=NaN;

%calc RA (s/m) based on EQ 2 on dan li paper (doi:10.1029/2018JG004401)
constant_cp=1000;%specific heat of air is 1 J·g-1·K-1
constant_rho_air=1.225;%kg/m³ 
data1_in_daily_RA=(data1_in_daily_TSK-data1_in_daily_T2)./(data1_in_daily_HFX./(constant_cp*constant_rho_air));%s/m

%LH (W/m2 or J/m2.s) (calc and use ETA for urban which reflect veg (pervious) portion of urban gridcells
[ data1_in_daily_LH           ,~,varunit_LH    ,vardesc_LH    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'LH'           );
data1_in_daily_LH =squeeze(nanmean(data1_in_daily_LH ,3));
[ data1_in_daily_LH_URB       ,~,~             ,~             ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'LH_URB'       );
data1_in_daily_LH_URB =squeeze(nanmean(data1_in_daily_LH_URB ,3));
data1_in_daily_ETA=(data1_in_daily_LH-data1_in_daily_LH_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data1_in_daily_LH_URB
data1_in_daily_LH(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data1_in_daily_ETA(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data1_in_daily_ETA
data1_in_daily_LH           (~CA_poly_on_in_withOcean)=NaN;

%Irrgation (for referece; not used in any calculations) 
[ data1_in_daily_NOAHRES    ,~,varunit_NOAHRES,vardesc_NOAHRES,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'NOAHRES'      );
data1_in_daily_NOAHRES       =squeeze(nanmean(data1_in_daily_NOAHRES     ,3));
data1_in_daily_NOAHRES           (~CA_poly_on_in_withOcean)=NaN;
data1_in_daily_NOAHRES=data1_in_daily_NOAHRES./(1-DATA_FRC_URB2D);%convert back to reflect veg (pervious) portion of urban gridcells (not gridcell level)
varunit_NOAHRES='mm/day';
Irr1=data1_in_daily_NOAHRES;clear data1_in_daily_NOAHRES

%surface available energy: LH + SH (W/m2) or (J/m2.s)
data1_in_daily_Rn=data1_in_daily_LH+data1_in_daily_HFX;

%Vapor pressure deficit (kPa)############note the unit is wrong to be Pa
[ data1_in_daily_VPD          ,~,varunit_VPD   ,vardesc_VPS   ,~,~,~] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'VPD'           );
data1_in_daily_VPD          =squeeze(nanmean(data1_in_daily_VPD          ,3));
data1_in_daily_VPD          (~CA_poly_on_in_withOcean)=NaN;
varunit_VPD='kPa';

%gradient of the saturation vapour pressure with respect to temperature (kPa/C)
[ data1_in_daily_GSVPT        ,~,varunit_GSVPT ,vardesc_GSVPT ,~,~,~] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'GSVPT'         );
data1_in_daily_GSVPT        =squeeze(nanmean(data1_in_daily_GSVPT        ,3));
data1_in_daily_GSVPT        (~CA_poly_on_in_withOcean)=NaN;
%% Calculate the partial differentials (based on  Hydrologic implications of vegetation response to elevated CO2 in climate projections
%                                       Yuting Yang" "1,2,6*, Michael L. Roderick" "1,2,3*, Shulei Zhang4, Tim R. McVicar2,5 and Randall J. Donohue2,5
%                                       https://doi.org/10.1038/s41558-018-0361-0
%constants
%psychrometric constant
gamma=0.0667;%(kPa/C)
%latent heat of vaporization : https://en.wikipedia.org/wiki/Latent_heat
lambda=2264.705/1000;%(MJ/kg) 

%variable for partial differentials calculation
%convert Rnstar and ET units to (mm/12hrDay) for potential ET calculation
ET1    =data1_in_daily_LH/(lambda*1000000)*(12*60*60);%(mm/12hrDay): J/m2.s / J/kg = kg/m2.s or mm/s x 12x60x60s/12hrDay = mm/12hrDay
RnStar1=data1_in_daily_Rn/(lambda*1000000)*(12*60*60);%(mm/12hrDay): same as above
Ta1    =data1_in_daily_T2   ;%(C)
s1     =data1_in_daily_GSVPT;%(kPa/C)
D1     =data1_in_daily_VPD  ;%(kPa)
ra1    =data1_in_daily_RA   ;%(s/m)
%Calculate rs by reversing Penman Monteith given in a working example that the author (Yuting) sent (a word doc)
rs1=(s1.*RnStar1.*ra1+(12*60*60).*((2.1458.*gamma)./(Ta1+273.15)).*D1)./(ET1.*gamma)-s1.*ra1./gamma-ra1;%(mm/12hrDay)
rs1(ET1==0)=NaN;%not to have inft 
rs1(DATA_LU_INDEX==17)=0;%Penman-OW: simplified Penman?Monteith over open-water used to calculate EP (rs = 0 and an aerodynamic roughness of the surface of 0.00137 m).
% alternative rs
% %calc RS (s/m) based on EQ 12 on dan li paper (doi:10.1029/2018JG004401)
% %calc qs at ts
% % saturation vapour pressure (kPa)
% data1_in_daily_es=0.6108*exp(17.27.*data1_in_daily_TSK./(data1_in_daily_TSK+237.3));clear data1_in_daily_TSK
% %load total pressure (Pa)?????????? is not at 2 m?
% [ data1_in_daily_PSFC         ,~,varunit_PSFC  ,vardesc_PSFC  ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'PSFC'         );
% data1_in_daily_PSFC =squeeze(nanmean(data1_in_daily_PSFC ,3));
% data1_in_daily_PSFC=data1_in_daily_PSFC/1000;%covert to kPa
% %calculate vapor pressure based on specific humidity: Shuttleworth book: q=0.622*e/P
% data1_in_daily_qs=0.622*(data1_in_daily_es./data1_in_daily_PSFC);clear data1_in_daily_es data1_in_daily_PSFC
% %LOAD atmospheric specific humidity (kg/kg)
% [ data1_in_daily_Q2           ,~,varunit_Q2    ,vardesc_Q2    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_1,'Q2'           );
% data1_in_daily_Q2 =squeeze(nanmean(data1_in_daily_Q2 ,3));
% %CALC RS based on EQ 12 on dan li paper (doi:10.1029/2018JG004401)
% lambda=2264.705*1000;%latent heat of vaporization (J/kg): https://en.wikipedia.org/wiki/Latent_heat
% constant_rho_air=1.225;%the air density (kg/m³) 
% data1_in_daily_RS=((constant_rho_air.*lambda).*(data1_in_daily_qs-data1_in_daily_Q2)./data1_in_daily_LH) - data1_in_daily_RA;%s/m
% rs1    =data1_in_daily_RS   ;%(s/m)
%% Calculate partial differentials
%convert Rnstar unit to MJ/m2.s
RnStar1_u=RnStar1/(12*60*60)*lambda;%(MJ/m2.s): mm/12hrDay /(12*60*60) = mm/s or kg/m2.s x MJ/kg = MJ/m2.s
%convert ET unit to kg/m2.s
ET1_u =ET1 /(12*60*60);%(kg/m2.s): mm/12hrDay / (12*60*60) = mm/s or kg/m2.s

%term 1: RHOa x Cp (MJ/C.m3): based on a working example that the author (Yuting) sent (a word doc)
term1=2.1458*gamma*lambda./(Ta1+273.15);%MJ/C.m3
%term2
term2=s1+gamma.*(1+rs1./ra1);%(kPa/C)+( (kPa/C) x () = (kPa/C)

% partial differentials
pd_E_RnStar=(s1./term2).*(1./lambda);%((kPa/C) / (kPa/C)) x 1/(MJ/kg) = kg/MJ
pd_E_D=term1./(ra1.*term2).*(1./lambda);%(MJ/C.m3) / (s/m x kPa/C) x (kg/MJ) = MJ/C.m3 / s.kPa/m.C x kg/MJ= (kg/m2.s.kPa)
pd_E_ra=(gamma.*rs1.*(s1.*RnStar1_u+term1.*D1./ra1))./(lambda.*ra1.^2.*term2.^2) - (term1.*D1)./(lambda.*ra1.^2.*term2);%kg/m.s2 - kg/m2.s.kPa x kPa / s/m = kg/m.s2
pd_E_rs=(-1*gamma.*(s1.*RnStar1_u+term1.*D1./ra1))./(lambda.*ra1.*term2.^2);% (kPa/C x (kPa.MJ/m2.s.C) / (MJ/kg x s/m x kPa2/C2) = kPa2.MJ/m2.s.C2 / MJ.s.kPa2/kg.m.C2 = (kg/m.s2)
pd_E_s=RnStar1_u./(lambda.*term2) - (s1.*RnStar1_u+term1.*D1./ra1)./(lambda.*term2.^2);% MJ/m2.s/(MJ/kg x kPa/C) + ... = kg.C/m2.s.kPa + ...

% Calculate relative sensitivity of E to change in Rn*, D, rs, ra and s, respectively (that is, relative changes in E induced by relative changes in the forcing variables)
pd_E_RnStar_ra=pd_E_RnStar.*RnStar1_u./ET1_u; %kg/MJ x  MJ/m2.s = kg/m2.s / kg/m2.s
pd_E_D_ra     =pd_E_D.*D1./ET1_u; %kg/m2.s.kPa x kPa = kg/m2.s / kg/m2.s
pd_E_ra_ra    =pd_E_ra.*ra1./ET1_u; %kg/m.s2 x s/m = kg/m2.s / kg/m2.s
pd_E_rs_ra    =pd_E_rs.*rs1./ET1_u; %kg/m.s2 x s/m = kg/m2.s / kg/m2.s
pd_E_s_ra     =pd_E_s.*s1./ET1_u; %kg/m2.s / kg/m2.s

% statistics for relative sensitivities
stat_pd_E_RnStar_ra = pouya_time_series2( pd_E_RnStar_ra*100, data_LU_INDEX, LULC, LULCl );
stat_pd_E_D_ra      = pouya_time_series2( pd_E_D_ra     *100, data_LU_INDEX, LULC, LULCl );
stat_pd_E_ra_ra     = pouya_time_series2( pd_E_ra_ra    *100, data_LU_INDEX, LULC, LULCl );
stat_pd_E_rs_ra     = pouya_time_series2( pd_E_rs_ra    *100, data_LU_INDEX, LULC, LULCl );
stat_pd_E_s_ra      = pouya_time_series2( pd_E_s_ra     *100, data_LU_INDEX, LULC, LULCl );
%%
for run=3%[4,2,5,3]
    eval(sprintf('runname_2=runname%d',run))
    eval(sprintf('plot_runname_2=plot_runname%d',run))

%% load future scenario data
%TSK (C) (calc and use ETA for urban which reflect veg (pervious) portion of urban gridcells
[ data2_in_daily_TSK          ,~,varunit_TSK   ,vardesc_TSK   ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'TSK'          );
data2_in_daily_TSK =squeeze(nanmean(data2_in_daily_TSK ,3));
[ data2_in_daily_TS_URB       ,~,~             ,~             ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'TS_URB'       );
data2_in_daily_TS_URB =squeeze(nanmean(data2_in_daily_TS_URB ,3));
data2_in_daily_T1=(data2_in_daily_TSK-data2_in_daily_TS_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data2_in_daily_TS_URB
data2_in_daily_TSK(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data2_in_daily_T1(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data2_in_daily_T1

%HFX (W/m2 or J/m2.s) (calc and use ETA for urban which reflect veg (pervious) portion of urban gridcells
[ data2_in_daily_HFX          ,~,varunit_HFX   ,vardesc_HFX   ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'HFX'          );
data2_in_daily_HFX =squeeze(nanmean(data2_in_daily_HFX ,3));
[ data2_in_daily_SH_URB       ,~,varunit_SH_URB,vardesc_SH_URB,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'SH_URB'       );
data2_in_daily_SH_URB =squeeze(nanmean(data2_in_daily_SH_URB ,3));
data2_in_daily_SHEAT=(data2_in_daily_HFX-data2_in_daily_SH_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data2_in_daily_SH_URB
data2_in_daily_HFX(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data2_in_daily_SHEAT(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data2_in_daily_SHEAT

%T2 (C)
[ data2_in_daily_T2           ,~,varunit_T2    ,vardesc_T2    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'T2'           );
data2_in_daily_T2 =squeeze(nanmean(data2_in_daily_T2 ,3));
data2_in_daily_T2           (~CA_poly_on_in_withOcean)=NaN;

%calc RA (s/m) based on EQ 2 on dan li paper (doi:10.1029/2018JG004401)
constant_cp=1000;%specific heat of air is 1 J·g-1·K-1
constant_rho_air=1.225;%kg/m³ 
data2_in_daily_RA=(data2_in_daily_TSK-data2_in_daily_T2)./(data2_in_daily_HFX./(constant_cp*constant_rho_air));%s/m

%LH (W/m2 or J/m2.s) (calc and use ETA for urban which reflect veg (pervious) portion of urban gridcells
[ data2_in_daily_LH           ,~,varunit_LH    ,vardesc_LH    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'LH'           );
data2_in_daily_LH =squeeze(nanmean(data2_in_daily_LH ,3));
[ data2_in_daily_LH_URB       ,~,~             ,~             ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'LH_URB'       );
data2_in_daily_LH_URB =squeeze(nanmean(data2_in_daily_LH_URB ,3));
data2_in_daily_ETA=(data2_in_daily_LH-data2_in_daily_LH_URB.*DATA_FRC_URB2D)./(1-DATA_FRC_URB2D);clear data2_in_daily_LH_URB
data2_in_daily_LH(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26)=data2_in_daily_ETA(DATA_LU_INDEX==24|DATA_LU_INDEX==25|DATA_LU_INDEX==26);clear data2_in_daily_ETA
data2_in_daily_LH           (~CA_poly_on_in_withOcean)=NaN;

%Irrgation (for referece; not used in any calculations) 
[ data2_in_daily_NOAHRES    ,~,varunit_NOAHRES,vardesc_NOAHRES,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'NOAHRES'      );
data2_in_daily_NOAHRES       =squeeze(nanmean(data2_in_daily_NOAHRES     ,3));
data2_in_daily_NOAHRES           (~CA_poly_on_in_withOcean)=NaN;
data2_in_daily_NOAHRES=data2_in_daily_NOAHRES./(1-DATA_FRC_URB2D);%convert back to reflect veg (pervious) portion of urban gridcells (not gridcell level)
varunit_NOAHRES='mm/day';
Irr2=data2_in_daily_NOAHRES;clear data2_in_daily_NOAHRES

%surface available energy: LH + SH (W/m2) or (J/m2.s)
data2_in_daily_Rn=data2_in_daily_LH+data2_in_daily_HFX;

%Vapor pressure deficit (kPa)############note the unit is wrong to be Pa
[ data2_in_daily_VPD          ,~,varunit_VPD   ,vardesc_VPS   ,~,~,~] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'VPD'           );
data2_in_daily_VPD          =squeeze(nanmean(data2_in_daily_VPD          ,3));
data2_in_daily_VPD          (~CA_poly_on_in_withOcean)=NaN;
varunit_VPD='kPa';

%gradient of the saturation vapour pressure with respect to temperature (kPa/C)
[ data2_in_daily_GSVPT        ,~,varunit_GSVPT ,vardesc_GSVPT ,~,~,~] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'GSVPT'         );
data2_in_daily_GSVPT        =squeeze(nanmean(data2_in_daily_GSVPT        ,3));
data2_in_daily_GSVPT        (~CA_poly_on_in_withOcean)=NaN;
%%
%constants
%psychrometric constant
gamma=0.0667;%(kPa/C)
%latent heat of vaporization : https://en.wikipedia.org/wiki/Latent_heat
lambda=2264.705/1000;%(MJ/kg) 
%variable for partial differentials calculation
%convert Rnstar and ET units to (mm/12hrDay) for potential ET calculation
ET2    =data2_in_daily_LH/(lambda*1000000)*(12*60*60);%(mm/12hrDay): J/m2.s / J/kg = kg/m2.s or mm/s x 12x60x60s/12hrDay = mm/12hrDay
RnStar2=data2_in_daily_Rn/(lambda*1000000)*(12*60*60);%(mm/12hrDay): same as above
D2     =data2_in_daily_VPD  ;%(kPa)
s2     =data2_in_daily_GSVPT;%(kPa/C)
ra2    =data2_in_daily_RA   ;%(s/m)
Ta2    =data2_in_daily_T2   ;%(C)
%Calculate rs (s/m) by reversing Penman Monteith given in a working example that the author (Yuting) sent (a word doc)
rs2=(s2.*RnStar2.*ra2+(12*60*60).*((2.1458.*gamma)./(Ta2+273.15)).*D2)./(ET2.*gamma)-s2.*ra2./gamma-ra2;%(s/m)
rs2(ET2==0)=NaN;%not to have inft 
rs2(DATA_LU_INDEX==17)=0;%Penman-OW: simplified Penman?Monteith over open-water used to calculate EP (rs = 0 and an aerodynamic roughness of the surface of 0.00137 m).
%alternative rs
% %calc RS (s/m) based on EQ 12 on dan li paper (doi:10.1029/2018JG004401)
% %calc qs at ts
% % saturation vapour pressure (kPa)
% data2_in_daily_es=0.6108*exp(17.27.*data2_in_daily_TSK./(data2_in_daily_TSK+237.3));clear data2_in_daily_TSK
% %load total pressure (Pa)?????????? is not at 2 m?
% [ data2_in_daily_PSFC         ,~,varunit_PSFC  ,vardesc_PSFC  ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'PSFC'         );
% data2_in_daily_PSFC =squeeze(nanmean(data2_in_daily_PSFC ,3));
% data2_in_daily_PSFC=data2_in_daily_PSFC/1000;%covert to kPa
% %calculate vapor pressure based on specific humidity: Shuttleworth book: q=0.622*e/P
% data2_in_daily_qs=0.622*(data2_in_daily_es./data2_in_daily_PSFC);clear data2_in_daily_es data2_in_daily_PSFC
% %LOAD atmospheric specific humidity (kg/kg)
% [ data2_in_daily_Q2           ,~,varunit_Q2    ,vardesc_Q2    ,~,~,~ ] = pouya_Processedwrfout_matlab_read_daily_daytime( foldername,runname_2,'Q2'           );
% data2_in_daily_Q2 =squeeze(nanmean(data2_in_daily_Q2 ,3));
% %CALC RS based on EQ 12 on dan li paper (doi:10.1029/2018JG004401)
% lambda=2264.705*1000;%latent heat of vaporization (J/kg): https://en.wikipedia.org/wiki/Latent_heat
% constant_rho_air=1.225;%the air density (kg/m³) 
% data2_in_daily_RS=((constant_rho_air.*lambda).*(data2_in_daily_qs-data2_in_daily_Q2)./data2_in_daily_LH) - data2_in_daily_RA;%s/m
% rs2    =data2_in_daily_RS   ;%(s/m)
%% Calculate attribution terms
%convert Rnstar unit to MJ/m2.s
RnStar2_u=RnStar2/(12*60*60)*lambda;%(MJ/m2.s): mm/12hrDay /(12*60*60) = mm/s or kg/m2.s x MJ/kg = MJ/m2.s

%atrribution terms
att_RnStar=pd_E_RnStar.*(RnStar2_u-RnStar1_u).*(12*60*60); %kg/MJ x MJ/m2.s = kg/m2.s or mm/s x(12*60*60) = mm/12hrDay 
att_D     =pd_E_D     .*(D2       -D1       ).*(12*60*60); %kg/m2.s.kPa x kPa = kg/m2.s or mm/s x(12*60*60) = mm/12hrDay
att_ra    =pd_E_ra    .*(ra2      -ra1      ).*(12*60*60); %kg/m.s2 x s/m = kg/m2.s or mm/s x(12*60*60) = mm/12hrDay
att_rs    =pd_E_rs    .*(rs2      -rs1      ).*(12*60*60); %kg/m.s2 x s/m = kg/m2.s or mm/s x(12*60*60) = mm/12hrDay
att_s     =pd_E_s     .*(s2       -s1       ).*(12*60*60); %kg.C/m2.s.kPa x kPa/C = kg/m2.s or mm/s x(12*60*60) = mm/12hrDay

% statistics for atrribution terms (mm/12hr daytime)
stat_att_RnStar = pouya_time_series2( att_RnStar, data_LU_INDEX, LULC, LULCl );
stat_att_D      = pouya_time_series2( att_D     , data_LU_INDEX, LULC, LULCl );
stat_att_ra     = pouya_time_series2( att_ra    , data_LU_INDEX, LULC, LULCl );
stat_att_rs     = pouya_time_series2( att_rs    , data_LU_INDEX, LULC, LULCl );
stat_att_s      = pouya_time_series2( att_s     , data_LU_INDEX, LULC, LULCl );

stat_tot_ET     = pouya_time_series2((ET2-ET1), data_LU_INDEX, LULC, LULCl );%kg/m2.s or mm/s x(12*60*60) = mm/12hrDay

%save results
savename=sprintf('%s/Data_partial_differentials_%s-%s_v11.mat',savefolder,plot_runname_2,plot_runname_1);
save(savename,'pd_E_*','stat_*','att_*','Irr1','ET1','RnStar1','D1','s1','ra1','Ta1','rs1','Irr2','ET2','RnStar2','D2','s2','ra2','Ta2','rs2','LULC*','data_*','savefolder','-v7.3')

end


