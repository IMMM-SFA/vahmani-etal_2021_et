_your zenodo badge here_

# vahmani-etal_2021_et

**Will anthropogenic warming increase Evapotranspiration? Examining Irrigation Water Demand Implications of Climate Change in California**

P. Vahmani<sup>\*1</sup>, A. D. Jones<sup>1</sup>, and D. Li<sup>2</sup>

<sup>1 </sup> Lawrence Berkeley National Laboratory, One Cyclotron Road, Berkeley, CA 94720, USA.

<sup>2 </sup> Department of Earth and Environment, Boston University, Boston, MA, USA.

\* corresponding author: pvahmani@lbl.gov

## Abstract
Climate modeling studies and observations do not fully agree on the implications of anthropogenic warming for evapotranspiration (ET), a major component of the water cycle and driver of irrigation water demand. Here we use California as a testbed to assess the ET impacts of changing atmospheric conditions induced by climate change on irrigated systems. Our analysis of irrigated agricultural and urban regions shows that warmer atmospheric temperatures have minimal implications for ET rates and irrigation water demands--about one percent change per degree Celsius warming (~1%°C<sup>-1</sup>). By explicitly modeling irrigation, we control for the confounding effect of climate-driven soil moisture changes and directly estimate water demand implications. Our attribution analysis of the drivers of ET response to global anthropogenic warming shows that as the atmospheric temperature and vapor pressure deficit depart from the ideal conditions for transpiration, regulation of stomata resistance by stressed vegetation almost completely offsets the expected increase in ET rates that would otherwise result from abiotic processes alone. We further show that anthropogenic warming of the atmosphere has minimal implications for mean relative humidity (<1.7%°C<sup>-1</sup>) and the surface available energy (<0.2%°C<sup>-1</sup>), which are critical drivers of ET. This study corroborates the growing evidence that plant physiological changes moderate the degree to which changes in potential ET are realized as actual ET.

## Journal reference
P. Vahmani P., A. D. Jones, and D. Li, 2021, Will anthropogenic warming increase Evapotranspiration? Examining Irrigation Water Demand Implications of Climate Change in California, Earth's Future.

## Data reference

### Input data
National Centers for Environmental Prediction/National Weather Service/NOAA/U.S. Department of Commerce (2005), NCEP North American Regional Reanalysis (NARR), https://rda.ucar.edu/datasets/ds608.0/, Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory, Boulder, Colo.

### Output data
WRF-UCM output data at 1.5 km resolution and supporting files can be accessed via [Globus endpoint](https://app.globus.org/file-manager?origin_id=86d8b02e-5948-11ec-b2c1-1b99bfd4976a&origin_path=%2F).

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| WRF | 3.6.1 | https://github.com/hanschen/WRFV3 | NA |

## Reproduce my experiment
1. Follow the instruction at https://www2.mmm.ucar.edu/wrf/OnLineTutorial/ to Download and compile WRFv3.6.1
2. Run WRF-UCM over CA for 15 years (2001-2015):
    - Use NARR reanalysis data (available at https://rda.ucar.edu/datasets/ds608.0/) as initial and boundary conditions. The workflow is available at https://www2.mmm.ucar.edu/wrf/OnLineTutorial/)
    - To reproduce the configuration used for this study use the corresponding namelist.input from the folder [1.namelist](1.namelist/) and the table of parameters from the folder [2.tables_edited](2.tables_edited/)
3. Repeat step 2 for three different scenarios as described in the manuscript. The WRF folders (one for each scenario) with the inputs, parameter tables, and namelists are available via [Globus endpoint](https://app.globus.org/file-manager?origin_id=86d8b02e-5948-11ec-b2c1-1b99bfd4976a&origin_path=%2F):
    - S31_base_ens01_WRFRun_2001.tar
    - S31_cnrm_ens01_WRFRun_2001.tar
    - S31_hadg_ens01_WRFRun_2001.tar
4. The attribution analysis (Figure 3) can be reproduced using the spreadsheet:
    - [All_BarPlot_Attribution_and_PrecChange_v2.xlsx](All_BarPlot_Attribution_and_PrecChange_v2.xlsx)

## Reproduce my figures
Use the scripts from the folder [4.Figures](4.Figures/) to reproduce the figures used in this publication:

| Script Name | Description | How to Run |
| --- | --- | --- |
| M01_Variable_Extraction_for_hpc_V2.m | Script to convert NetCDF files to MATLAB files | `matlab -nodisplay -nosplash -nodesktop -r "run('M01_Variable_Extraction_for_hpc_V2.m');exit;"` |
| aM01_wrfout_plot_spatial_v4.m | Script to create spatial maps in Figures 2, 4, S1, S2, S4 | `matlab -nodisplay -nosplash -nodesktop -r "run(‘aM01_wrfout_plot_spatial_v4.m');exit;"` |
| aM0_wrfout_attribution_of_changes_in_ET_v11.m | Script to create intermediate data on partial differentials presented in the attribution analysis | `matlab -nodisplay -nosplash -nodesktop -r "run('aM0_wrfout_attribution_of_changes_in_ET_v11.m');exit;"` |
| aM1_change_spatial.m<br><br><br>aM2_change_bar.m<br><br><br>aM3_attribution_bar.m | Scripts to use intermediate data from above to create Figures 3 and S3 | `matlab -nodisplay -nosplash -nodesktop -r "run(aM1_change_spatial.m');exit;"`<br><br>`matlab -nodisplay -nosplash -nodesktop -r "run(aM1_change_bar.m');exit;"`<br><br>`matlab -nodisplay -nosplash -nodesktop -r "run(aM1_attricution_bar.m');exit;"` |

