######## *** FUNCTION TO BUILD PARAMETER LISTS FOR GIVEN COUNTRY  *** ##########
# Builds a list of input parameters for each country for the decision tree. 
# Selects the appropriate parameter samples from param.list based on country/region 
# and combines this with country specific parameter values from data.countries

make_params <- function(country,dl,n_samp){

  country_params <- list()
  
  #-------- select samples for dl$parameters based on country/region --------
  
  for (name in names(dl$samples)){                
    record <- dl$parameters[full_name==name,]  
    reg_class <- record[["region_class"]]   #region def. to map this parameter

    # if parameter common to all countries then just add samples to the list
    if (reg_class == ""){                        
      country_params[[name]] <- dl$samples[[name]]
    }
    
    # otherwise look up the appropriate regional parameter value using region_class
    else if (record[["region"]] == dl$regions.dt[iso3==country, get(reg_class)]){  
      
      # remove the region suffix from the parameter name to give the generic name used in the tree
      generic_name <- gsub("\\.\\w*$","",name)
      country_params[[generic_name]] <- dl$samples[[name]]
    }
  }
  
  #-------- now add other country specific parameters --------
  
  #-------- GDP per capita --------
  country_params[["gdp_per_cap"]] <- dl$countries.dt[iso3==country, wb_gdp_per_capita_usd_2020]
  
  #-------- vaccine coverage --------
  country_params[["vac_cov_assumption"]] <- "ANC4" # default vaccine coverage assumption 
  country_params[["p.anc4"]] <- dl$countries.dt[iso3==country, p.anc4.ihme.2019]
  country_params[["p.anc1"]] <- dl$countries.dt[iso3==country, p.anc1.ihme.2019]
  country_params[["p.preterm_covered"]] <- dl$countries.dt[iso3==country, p.preterm_covered]
  
  #-------- life expectancy --------
  country_params[["life_expectancy"]] <- dl$countries[iso3==country,][["life_exp_birth_2020"]]
  
  #-------- maternal colonisation prevalence --------
  country_params[["p.mat_col"]] <- get_country_data(dl$mat_col.dt,country,dl$regions.dt,n_samp)
  
  #-------- stillbirth --------
  # overall risk of stillbirth and livebirth
  risk.stillbirth <- dl$countries.dt[iso3==country, Stillbirths_2020] / dl$countries.dt[iso3==country, Births_2020]
  
  country_params[["p.livebirth"]] <- 1 - risk.stillbirth 
  
  # proportion of stillbirth attributed to gbs
  prop.gbs_stillbirth <- get_country_data(dl$prop_gbs_stillbirth.dt,country,dl$regions.dt,n_samp)
  
  country_params[["risk.stillbirth.no_col"]] <- risk.stillbirth * (1 - prop.gbs_stillbirth)
  country_params[["risk.stillbirth.col"]] <-   country_params[["risk.stillbirth.no_col"]] + (risk.stillbirth * prop.gbs_stillbirth) / country_params[["p.mat_col"]]
 
  #-------- preterm birth -------- 
  # call or2r which returns adjusted risk among exposed / unexposed 
  adjusted_preterm <- or2r(
    dl$countries.dt[iso3==country, Preterm_2014],                      # preterm risk among all live births
    country_params[["p.mat_col"]],                                     # mat. col. among live births
    get_country_data(dl$OR_preterm.dt,country,dl$regions.dt,n_samp)    # OR for preterm given mat. col.
  )
  
  country_params[["risk.preterm.no_col"]] <- adjusted_preterm$r0
  country_params[["risk.preterm.col"]] <- adjusted_preterm$re
  
  #-------- risk of eod --------
  country_params[["p.eod.col"]] <- get_country_data(dl$risk_gbs.dt ,country,dl$regions.dt,n_samp)
  
  #-------- risk of lod --------
  # proportion of gbs that lod
  country_params[["prop.lod"]]   <- get_country_data(dl$prop_lod.dt ,country,dl$regions.dt,n_samp)
  
  #-------- case fatality risk --------
  
  country_params[["coverage_SBA"]] <- dl$countries.dt[iso3==country, coverage_SBA]/100
  country_params[["p.gbs_death.eod"]] <- get_country_data(dl$cfr_eod.dt ,country,dl$regions.dt,n_samp)
  country_params[["p.gbs_death.lod"]] <- get_country_data(dl$cfr_lod.dt ,country,dl$regions.dt,n_samp)
  
  #-------- meningitis --------
  country_params[["p.men.eod"]] <- get_country_data(dl$prop_eod_men.dt ,country,dl$regions.dt,n_samp)
  
  country_params[["p.men.lod"]] <- get_country_data(dl$prop_lod_men.dt ,country,dl$regions.dt,n_samp)
  
  #-------- NDI --------
  # after meningitis
  p.ndi.men <- get_country_data(dl$risk_ndi_men.dt ,country,dl$regions.dt,n_samp)
  p.ndi_modsev.men <- get_country_data(dl$risk_ndi_men_modsev.dt ,country,dl$regions.dt,n_samp)

  # These are now adjusted to give GBS attributable NDI
  country_params[["p.mild_ndi.men"]] <- p.ndi.men -   p.ndi_modsev.men - country_params[["p.mild_ndi_baseline"]]
  country_params[["p.mod_ndi.men"]]  <- (p.ndi_modsev.men - country_params[["p.mod_sev_ndi_baseline"]]) * (1 - country_params[["p.mod_sev_ndi_sev"]])
  country_params[["p.sev_ndi.men"]]  <- (p.ndi_modsev.men - country_params[["p.mod_sev_ndi_baseline"]]) * country_params[["p.mod_sev_ndi_sev"]]
  country_params[["p.no_ndi.men"]] <- 1 - country_params[["p.mild_ndi.men"]] - country_params[["p.mod_ndi.men"]] - country_params[["p.sev_ndi.men"]] 
  
  # absolute any and mod/sev ndi risk after meningitis
  country_params[["p.modsev_ndi_abs.men"]] <- p.ndi_modsev.men
  country_params[["p.any_ndi_abs.men"]]    <- p.ndi.men
  
  # after sepsis
  p.ndi.sep <- get_country_data(dl$risk_ndi_sep.dt ,country,dl$regions.dt,n_samp)
  p.ndi_modsev.sep <- get_country_data(dl$risk_ndi_sep_modsev.dt ,country,dl$regions.dt,n_samp)
  
  # These are now adjusted to give GBS attributable NDI
  country_params[["p.mild_ndi.sep"]] <- p.ndi.sep - p.ndi_modsev.sep - country_params[["p.mild_ndi_baseline"]]
  country_params[["p.mod_ndi.sep"]]  <- (p.ndi_modsev.sep - country_params[["p.mod_sev_ndi_baseline"]]) * (1 - country_params[["p.mod_sev_ndi_sev"]])
  country_params[["p.sev_ndi.sep"]]  <- (p.ndi_modsev.sep - country_params[["p.mod_sev_ndi_baseline"]]) * country_params[["p.mod_sev_ndi_sev"]]
  country_params[["p.no_ndi.sep"]] <- 1 - country_params[["p.mild_ndi.sep"]] - country_params[["p.mod_ndi.sep"]] - country_params[["p.sev_ndi.sep"]] 
  
  # absolute any and mod/sev ndi risk after sepsis
  country_params[["p.modsev_ndi_abs.sep"]] <- p.ndi_modsev.sep
  country_params[["p.any_ndi_abs.sep"]]    <- p.ndi.sep
  
  #-------- QALYs --------
  
  
  #-------- Costs --------
  
  # Acute costs 
  # log normal dist. based on meanlog and selog from fit to acute cost data
  
  country_params[["c.acute_hospital"]] <- rlnorm(
    n_samp,
    meanlog = dl$countries.dt[iso3==country, cost_acute_mn_log],
    sdlog = dl$countries.dt[iso3==country, cost_acute_se_log]
  )
  
  # Vax delivery costs 
  # log normal dist. based on meanlog and selog from fit to vax delivery data
  
  country_params[["c.vac_delivery"]] <- rlnorm(
    n_samp,
    meanlog = dl$countries.dt[iso3==country, cost_vax_del_mn_log],
    sdlog = dl$countries.dt[iso3==country, cost_vax_del_se_log]
  )
  
  return(country_params)
}

# function to extract country-specific data

get_country_data <- function(dt,country,regions.dt,n_samp){
  
  # if data is at country level
  if ("iso3"%in%names(dt)){ 
    r <- dt[iso3==country & sample<=n_samp, parameter_value]
  } 
  # if data is at regional level
  else if ("region"%in%names(dt)){ 
    reg_class <- dt[1,region_classification]
    reg <- regions.dt[iso3==country, get(reg_class)]
    r <- dt[region==reg & sample<=n_samp, parameter_value]
  } 
  # else assume a global sample
  else {r <- dt[sample<=n_samp, parameter_value]}
  
  # warning if n_samp is longer than number of samples available
  if (n_samp < length(r)){print("Warning no samples less than n_samp.")}
  
  return(r)
}
