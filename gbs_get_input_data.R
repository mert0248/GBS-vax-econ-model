######## *** FUNCTION TO READ IN INPUT DATA  *** ##########

source("gbs_helper_funcs.R")
require("data.table")

get_input_data <- function(path.io,
                           stochastic, 
                           n_samp, 
                           n_dist){
  
  # make list of data tables
  dl <- list(
    
    # placeholder for sampled data from data.parameters.csv
    samples = list(),
    # table of global and regional parameter distributions
    parameters = fread(paste0(path.io,"data.parameters.csv"),  
                       colClasses = c(mat_col = "character",
                                      sequelae = "character",
                                      prematurity = "character",
                                      onset = "character",
                                      iap_status = "character",
                                      region_class = "character",
                                      region = "character")),
    # table of country specific data
    countries.dt = as.data.table(
      fread(paste0(path.io,"data.countries.csv"))
    ),
    
    # table mapping country to region codes
    regions.dt = as.data.table(
      fread(paste0(path.io,"data.region.definitions.csv"))
    ),
    
    # maternal colonisation - from burden model
    mat_col.dt = as.data.table(
      readRDS(paste0(path.io,"GBS_colonisation_parameter_updated_18032021.rds"))
    ),

    # risk EOGBS if colonised - from burden model
    risk_gbs.dt = as.data.table(
      readRDS(paste0(path.io,"GBS_risk_givencolonisation_updated_18032021.rds"))
    ),
    
    # proportion of GBS that is LOGBS - from burden model
    prop_lod.dt = as.data.table(
      readRDS(paste0(path.io,"LOGBS_proportion_14022021.rds"))
    ),
    
    # case fatality risk for EOGBS - from burden model
    cfr_eod.dt = as.data.table(
      readRDS(paste0(path.io,"EO_CFR_14022021.rds"))
    ),
    
    # case fatality risk for LOGBS - from burden model
    cfr_lod.dt = as.data.table(
      readRDS(paste0(path.io,"LO_CFR_14022021.rds"))
    ),
    
    # proportion EOGBS that's meningitis - from burden model
    prop_eod_men.dt = as.data.table(
      readRDS(paste0(path.io,"overall_EO_meningitis_prop_04052021.rds"))
    ),
    
    # proportion LOGBS that's meningitis - from burden model
    prop_lod_men.dt = as.data.table(
      readRDS(paste0(path.io,"overall_LO_meningitis_prop_12022021.rds"))
    ),
    
    # risk of NDI following GBS meningitis - from burden model
    risk_ndi_men.dt = as.data.table(
      readRDS(paste0(path.io,"overall_men_any_NDI_07102021.rds"))
    ),
    
    # risk of mod/sev NDI following GBS meningitis - from burden model
    risk_ndi_men_modsev.dt = as.data.table(
      readRDS(paste0(path.io,"overall_men_modsev_NDI_07102021.rds"))
    ),
    
    # risk of NDI following GBS sepsis - from burden model
    risk_ndi_sep.dt = as.data.table(
      readRDS(paste0(path.io,"overall_sepsis_any_NDI_07102021.rds"))
    ),
    
    # risk of mod/sev NDI following GBS sepsis - from burden model
    risk_ndi_sep_modsev.dt = as.data.table(
      readRDS(paste0(path.io,"overall_sepsis_modsev_NDI_07102021.rds"))
    ),
    
    # proportion of all stillbirths due to GBS - from burden model
    prop_gbs_stillbirth.dt = as.data.table(
      readRDS(paste0(path.io,"stillbirth_GBS_proportion_10042021.rds"))
    ),
    
    # OR for preterm birth if colonised - from burden model
    OR_preterm.dt = as.data.table(
      readRDS(paste0(path.io,"OR_cohort_cross_sectional_preterm_15032021.rds"))
    )
    
  )
  
  # get samples from data.parameters, either fixed values or stochastic draws using
  # latin hyper-cube sampling
  dl[["samples"]] = get_samples(dl$parameters, stochastic, n_samp, n_dist)
  
  return(dl)
}