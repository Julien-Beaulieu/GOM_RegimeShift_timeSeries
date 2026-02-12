# GOM_RegimeShift_timeSeries
This repository is for all the codes and data used for time series and delayed synchrony identification. This includes all of the data cleaning steps.


## List of Script

#### 1st_2nd_derivative_NOAA_predators_UV_TS.R
This script is to calculate the first derivative of the time series
of subtidal predators directly feeding on intertidal preys. These 
species were identified in dataPrep_NOAA_intertidal_Web.R and 
NOAA-intertidal_interaction_web.R. The timeseries were fitted in 
fit_NOAA_UV_TS_linkers.R.
contributor: Julien Beaulieu

#### cross-correlation_inter_sub.R

This script is to calculate the cross-correlation between the temporal trends
of subtidal predators identified to feed on the intertidal species and each of their
subtidal preys. These were identified using binary connectance web with the GLOBI
interaction data base from the scripts dataPrep_NOAA-intertidal_Web.R and
NOAA-intertidal_interaction_web.R. Temporal trends in the intertidal are fitted
in time-sries_JB.R. Temporal trends in the subtidal ar fitted in fit_NOAA_UV_TS_linkers.R.
Contributor: Julien Beaulieu


#### dataPrep_NOAA-intertidal_Web.R

This script is to pull the interactions of the NOAA subtidal sp. 
with the intertidal species to identify key linkers. The binary connectance web
using this data is made in the script NOAA-intertidal_interaction_web.R.
This is the first step to propos a mechanism for RS propagation.
Contributor: Julien Beaulieu


#### derivative_calculation.R

This script is to calculate the FIRST derivative using forward difference.
This script is to calculate the FIRST derivative to know if there is a change
in intertidal communities when there is a regime shift. The script also produces
a table with the 1st deriv values during subtidal regime shifts and produce 1st
deriv figures.
Contributor: Julien Beaulieu

#### first_deriv_function.R

This script is a function to calculate First derivative with forward
difference the input is a dataframe with the year as first colums and the
draws posterior prediction as column. This is also the output of the extract 
function. 
contributor: Julien Beaulieu


#### fit_multivar_TS.R

This script is to fit multivariate time-series (GAM) on intertidal data and also 
subtidal data from NOAA and compare the rates of changes.
Contributor: Julien Beaulieu


#### fit_NOAA_UV_TS_linkers.R

This script is to fit a GAM on each subtidal specie from the NOAA
bottom trawl survey that has a direct interaction with intertidal taxa.
These species are identified using binary connectance web with the GLOBI
interaction data base from the scripts dataPrep_NOAA-intertidal_Web.R and
NOAA-intertidal_interaction_web.R
Contributor: Julien Beaulieu

#### format_sp_name_for_globi.R

This script is to make a function that change sp. name format
from genus_spcie to Genus specie. 


#### function_genusSp_2_GSp.R

This script makes a function to shift from genus specie to G. specie
Contributor: Julien Beaulieu

#### function_get_deriv_value.R

This script is to make a function that takes a first or second derivative 
and pulls out the derivatives values during regime shifts. The input is the 
output of the Deriv function with columns as: year, y_deriv_estim, sd, upper_CI,
and lower_CI.
Contributor: Julien Beaulieu


#### model_selection_mvTS.R

This script is to assess convergence and compare multivariate GAMs

#### data_prep_multivariate_time_series.R

This script is to prepare the data and ordinations to fit multivariate time-series 
(GAM) on intertidal data and also subtidal data from NOAA and compare the rates of changes.
Contributor: Julien beaulieu

#### noaa_AppIsland_sub_comparison

This script is to compare the species present in the NOAA bottom-trawl surveys 
to those observed in the Shouls subtidal survey (keen data)
Author: Julien Beaulieu


#### plot_deriv.R

This scrpt is to make a function to plot derivative trends. The input of this 
function is a data frame given by the Deriv function
Contributor: Julien Beaulieu


#### plot_heatmap_UV_TS.R

This script is to plot a heat map with the synchrony between changes in 
intertidal species abundances and subtidal regime shifts. The synchrony was
determined based on the codels and graphs in the script plot_UV_TS.R
Contributor: Julien Beaulieu


#### plot_MV_TS.R
This script si to plot time series on ordination axis.
The models used are those selected in model_selection_mvTS.R
the 2nd derivative values are from
contributor: Julien Beaulieu

#### plot-ordinations-JB.R

This script is to plot morinations. The score or the PCoA are from the script 
multivariate_time_series. The PCoA scores (input) are produced by the script 
multivariate_time_series
Contributor: Julien Beaulieu

#### plot-UV-TS.R

This script is to code graph for species-specific time series
The synchrony was determined by the script second_derivative_central_diff.R
contributor: Julien Beaulieu

#### pp_extract_function.R

This script is a funcition to extract spaghetti posterior prediction from 
brms models. It is used mainly in the script second_derivative_central-diff.R
Contributor: Julien Beaulieu

#### regime-shift-time_ID.R

This script is to calculate the year where each regime shift happenes (deriv max)
with a 95% CI.
Contributor: Julien Beaulieu

#### second_deriv_function.R

This script is a function to calculate the second derivative with central
difference the input is a dataframe with the year as first colums and the
draws posterior prediction as column. This is also the output of the extract 
function. 
contributor: Julien Beaulieu

#### second_derivative_central-diff.

This script is to calculate the second derivative using central difference.
This script is to calculate the second derivative to know if there is a change
in the speed of the change in intertidal communities (elbow) when there is 
a regime shift. The lscript also produces a table with the 2nd deriv values
during subtidal regime shifts and produce 2nd deriv figures.
Contributor: Julien Beaulieu

#### time_line_GOM_nat-History.R

This script is to produce a time line with the principal perturbation events
in the GOM from 1980 to 2020. This script also makes conceptual figures of the 
derivative.
Contributor: Julien Beaulieu

#### time-series.R

This script is to fit time series on the principal organism from the SEMs.
Contributor: Julien Beaulieu

#### time-series-model_selection.R

This script is to do model selection on the principal organism from the SEMs
Contributor: Julien Beaulieu

#### SEM_dataset_assembly_GCmeans.R

dataset assembly and cleaning.
Contributor: Nicole Knight


## List of data files

#### direct_interaction_int-sub.csv

List of interactions between subtidal species in the NOAA data and the key intertidal
taxa selected exclusively.

#### full_dataset_for_SEMs_groupcentered.csv

Intertidal data used as input for model fitting and ordinations

#### indirect_direct_interaction_int-sub.csv

Interactions between intertidal ans subtidal species including indirect interactions 
(interactions between species of the same zone indirectly connected to another zone).

#### KEEN_data.csv

Near subtidal data used to identify taxa present in the NOAA data that are commonly 
present in the higher subtidal

#### NOAA_trawl_data.csv

NOAA bottom trawl data used as input for models and ordinationd.

#### NOAA_trawl_locations.csv

Location of all six NOAA trawl transects around Appledore island that were selected.








