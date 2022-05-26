
# Decision tree model for the economic evaluation of GBS maternal vaccination

#########################   ***    SETTINGS   ***    ###########################

# set.seed = 1234
set.seed   = 6789
stochastic = TRUE # stochastic sampling from parameter distributions
n_samp     = 4000     # number of samples
n_dist     = 300      # number of parameter distributions used for LHS 
reg_type   = "sdg"    # region grouping for output "sdg" or "wb"

path.io <- "./output/" # path for outputs
path.dat <- "./data/"  # path for input data

#################   ***    REQUIRED PACKAGES   ***       #######################

source("gbs_helper_funcs.R")   # helper functions
source("gbs_get_input_data.R") # read in data files
source("gbs_make_params.R")    # generates country-specific parameters lists
source("gbs_eval_tree.R")      # evaluates tree

require("data.table") 
require("tidyverse")
require("ggplot2")
require("scales")
require("gridExtra")
require("ggsci")

#################  ***    LOOP OVER REGION GROUPINGS    ***   ##################

start_time <- Sys.time()
for (reg_type in c("sdg","wb")) { 

########################  ***   GET INPUT DATA   ***   #########################

dl <- get_input_data(path.dat, stochastic, n_samp, n_dist)

###############    ***   LIST OF COUNTRIES AND SCENARIOS  ***     ##############

# list of countries (iso3 codes) and scenarios
countries <- dl$countries[,iso3]
n_countries <- length(countries)

# load scenarios
scen_pars <- as.data.table(fread(paste0(path.dat,"data.scenarios.csv")))
n_scen <- max(scen_pars[,scen_id])
scen_pars.dt <- melt(
  scen_pars,
  id.vars=c("scen_id","base_id","vac_id","base_name","vac_name","scen_name"), 
  variable.name="par",
  value.name="val",
  variable.factor = FALSE
)

######################      ***   RUN MODEL   ***     ##########################

# list for results from countries x scenarios
country_results <- vector("list", n_scen * n_countries)

##### LOOP OVER COUNTRIES #####    
k = 0 # counter
for (i in 1:n_countries){
  country <- countries[i]
  cat("\r",i," of ", n_countries)

  #build baseline parameter lists for each country
  country_params <- make_params(country,dl,n_samp)
  pregnancies <- dl$countries.dt[iso3==country,Births_2020] # no. of pregnancies
  
  ##### LOOP OVER SCENARIOS #####
  for (s in 1:n_scen){
    k = k + 1 # increment counter
    params <- country_params                # reset parameters to baseline
    scen <- scen_pars.dt[scen_id==s]        # get changes for current scenario
    for (j in 1:length(scen[[1]])){         
      if(!is.na(scen[j,val])){              # don't update values that are NA
        params[[scen[j,par]]]=scen[j,val]   # update value in params using scen
      }
    }

    # evaluate model
    country_results[[k]] <- eval_tree(params)
    cols <- setdiff(names(country_results[[k]]),c("samp_id"))
    
    # multiply results by number of pregnancies
    country_results[[k]] <- country_results[[k]][
      ,lapply(.SD,
              function(x){x*pregnancies}),
      by=.(samp_id),
      .SDcols=cols
    ]
    
    # add scenario id, country details
    if (reg_type=="sdg"){
    region <- dl$regions.dt[iso3==country,sdg]
    } else if(reg_type =="wb") {
      region <- dl$regions.dt[iso3==country,wb_income]
    } else {print("Error: reg_type undefined")}
    base_id <- scen_pars[scen_id==s,base_id]
    vac_id <- scen_pars[scen_id==s,vac_id]
    country_results[[k]][, 
      c("iso3","scen_id","base_id","vac_id","region") := 
        .(country, s, base_id, vac_id, region)
    ]
    
  } ##### END OF SCENARIOS #####
} ##### END OF COUNTRY LOOP #####

# combine results into single data table
results.dt <- rbindlist(country_results)
country_results <- NULL # free memory

# melt
cols <- setdiff(names(results.dt),c("samp_id","scen_id","base_id","vac_id","iso3","region"))
results.dt <- melt(
  results.dt,
  id.vars=c("iso3","region","scen_id","base_id","vac_id","samp_id"), 
  measure.vars=cols,
  variable.name="name",
  value.name="val"
)

# calculate incremental values for vaccine scenarios vs corresponding no-vaccine
# baseline
print("Calculating incremental outcomes")

base.dt      <- results.dt[vac_id == 0] # no vaccine counter-factual scenarios
scenarios.dt <- results.dt[vac_id > 0] # vaccine scenarios
results.dt <- NULL # free memory

scenarios.dt <- scenarios.dt[base.dt, on=.(iso3,base_id,samp_id,name)][
  ,.(iso3, region, scen_id, base_id, vac_id, samp_id, name, val, inc=val-i.val)
]

##################         *** SUMMARISE OUTCOMES ***        ###################

# summarise for baseline (non-vaccination) scenarios
print("baseline summaries")

# global
base.global <- base.dt[,.(val=sum(val)), by=.(samp_id,scen_id,base_id,vac_id,name)][
  , .(val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975))),
  by=.(scen_id,base_id,vac_id,name)
]
base.global <- scen_pars[,.(scen_id,base_name,vac_name)][base.global, on=.(scen_id)] # add scenario names

# by region
base.regions <- base.dt[,.(val=sum(val)), by=.(samp_id,scen_id,base_id,vac_id,name,region)][
  , .(val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975))),
  by=.(scen_id,base_id,vac_id,name,region)
]
base.regions <- scen_pars[,.(scen_id,base_name,vac_name)][base.regions, on=.(scen_id)] # add scenario names

# by country
base.countries <- base.dt[,.(val=sum(val)), by=.(samp_id,scen_id,base_id,vac_id,name,iso3)][
  , .(val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975))),
  by=.(scen_id,base_id,vac_id,name,iso3)
]
base.countries <- scen_pars[,.(scen_id,base_name,vac_name)][base.countries, on=.(scen_id)] # add scenario names

# scenario summaries
print("scenario summaries")

# global
scenarios.global <- scenarios.dt[,.(val=sum(val), inc=sum(inc)), by=.(samp_id,scen_id,base_id,vac_id,name)][
  , .(
    val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975)),
    inc.mn=mean(inc), inc.md=median(inc), inc.lo95=quantile(inc,c(0.025)), inc.hi95=quantile(inc,c(0.975))
  ),
  by=.(scen_id,base_id,vac_id,name)
]
scenarios.global <- scen_pars[,.(scen_id,base_name,vac_name)][scenarios.global, on=.(scen_id)] # add scenario names

# by region
scenarios.regions <- scenarios.dt[,.(val=sum(val), inc=sum(inc)), by=.(samp_id,scen_id,base_id,vac_id,name,region)][
  , .(
    val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975)),
    inc.mn=mean(inc), inc.md=median(inc), inc.lo95=quantile(inc,c(0.025)), inc.hi95=quantile(inc,c(0.975))
  ),
  by=.(scen_id,base_id,vac_id,name,region)
]
scenarios.regions <- scen_pars[,.(scen_id,base_name,vac_name)][scenarios.regions, on=.(scen_id)] # add scenario names

# by country
scenarios.countries <- scenarios.dt[,.(val=sum(val), inc=sum(inc)), by=.(samp_id,scen_id,base_id,vac_id,name,iso3)][
  , .(
    val.mn=mean(val), val.md=median(val), val.lo95=quantile(val,c(0.025)), val.hi95=quantile(val,c(0.975)),
    inc.mn=mean(inc), inc.md=median(inc), inc.lo95=quantile(inc,c(0.025)), inc.hi95=quantile(inc,c(0.975))
  ),
  by=.(scen_id,base_id,vac_id,name,iso3)
]
scenarios.countries <- scen_pars[,.(scen_id,base_name,vac_name)][scenarios.countries, on=.(scen_id)] # add scenario names

base.dt <- NULL # free memory

##############         *** ICERS, NMB, THRESHOLD PRICE ***        ##############

print("icer, nmb, threshold price summaries")

icers.dt <- dcast(
  scenarios.dt[name %in% c(
    "qaly_loss.gbs_disease",
    "qaly_loss.gbs_preterm",
    "qaly_loss.gbs_stillbirth",
    "doses",
    "cost.vac_price",
    "cost.total")],
  iso3 + region + scen_id + samp_id ~ name,
  value.var = "inc"
)

# add country specific CET based on 0.5 and 1 x GDP per capita, and Woods/Ochalek thresholds
# drop countries for which there is no GDP per capita data

cets.dt <- rbind(
  dl$countries.dt[!is.na(wb_gdp_per_capita_usd_2020),.(iso3, cet = wb_gdp_per_capita_usd_2020, cet_name="1xGDPPC")],
  dl$countries.dt[!is.na(wb_gdp_per_capita_usd_2020),.(iso3, cet = cet_woods_ochalek, cet_name="Empirical")]
)

icers.dt <- merge(icers.dt,cets.dt,by="iso3",allow.cartesian = T)
cets.dt <- NULL

# create scenarios with and without valuation of stillbirth QALYs 
# (saves having to run these as separate scenarios in the main model)
sb.dt <- data.table(merge_id=1, with_SB=c("No","Yes"))
icers.dt <- merge(icers.dt[,merge_id:=1],sb.dt,by="merge_id",allow.cartesian = T)
icers.dt[,merge_id:=NULL]
sb.dt <- NULL

icers.dt[with_SB=="No", qaly_loss.gbs_stillbirth := 0]

icers.dt[, nmb := (-1 * cet * (qaly_loss.gbs_disease + qaly_loss.gbs_preterm + qaly_loss.gbs_stillbirth)) - cost.total]

# global
icers.global <- icers.dt[
  ,.(
    qaly_loss.gbs_disease = sum(qaly_loss.gbs_disease),
    qaly_loss.gbs_preterm = sum(qaly_loss.gbs_preterm),
    qaly_loss.gbs_stillbirth = sum(qaly_loss.gbs_stillbirth),
    cost.total = sum(cost.total),
    nmb = sum(nmb)
  ), 
  by=.(scen_id,with_SB,cet_name,samp_id)
][
  ,.(
    scen_id,with_SB,cet_name,samp_id, nmb,
    icer = cost.total / (-1 * (qaly_loss.gbs_disease + qaly_loss.gbs_preterm + qaly_loss.gbs_stillbirth))
  )
][
  ,.(
    nmb.mn=mean(nmb), 
    nmb.md=median(nmb),
    nmb.lo95=quantile(nmb,c(0.025)), 
    nmb.hi95=quantile(nmb,c(0.975)),
    icer.mn=mean(icer), 
    icer.md=median(icer),
    icer.lo95=quantile(icer,c(0.025)), 
    icer.hi95=quantile(icer,c(0.975))
  ),
  by=.(scen_id,with_SB,cet_name)
]
icers.global <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][icers.global, on=.(scen_id)] # add scenario names

# by region
icers.regions <- icers.dt[
  ,.(
    qaly_loss.gbs_disease = sum(qaly_loss.gbs_disease),
    qaly_loss.gbs_preterm = sum(qaly_loss.gbs_preterm),
    qaly_loss.gbs_stillbirth = sum(qaly_loss.gbs_stillbirth), 
    cost.total = sum(cost.total),
    nmb = sum(nmb)
  ), 
  by=.(region,scen_id,with_SB,cet_name,samp_id)
][
  ,.(
    region,scen_id,with_SB,cet_name,samp_id, nmb,
    icer = cost.total / (-1 * (qaly_loss.gbs_disease + qaly_loss.gbs_preterm + qaly_loss.gbs_stillbirth))
  )
][
  ,.(
    nmb.mn=mean(nmb), 
    nmb.md=median(nmb),
    nmb.lo95=quantile(nmb,c(0.025)), 
    nmb.hi95=quantile(nmb,c(0.975)),
    icer.mn=mean(icer), 
    icer.md=median(icer),
    icer.lo95=quantile(icer,c(0.025)), 
    icer.hi95=quantile(icer,c(0.975))
  ),
  by=.(region,scen_id,with_SB,cet_name)
]
icers.regions <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][icers.regions, on=.(scen_id)] # add scenario names

# by country - including vaccine threshold price

icers.dt[,vac_tp := (nmb + cost.vac_price) / doses]

icers.countries <- icers.dt[
  ,.(
    iso3, region, scen_id, with_SB, cet_name, nmb, vac_tp, samp_id,
    icer = cost.total / (-1 * (qaly_loss.gbs_disease + qaly_loss.gbs_preterm + qaly_loss.gbs_stillbirth))
  )
][
  ,.(
    nmb.mn=mean(nmb), 
    nmb.md=median(nmb),
    nmb.lo95=quantile(nmb,c(0.025)), 
    nmb.hi95=quantile(nmb,c(0.975)),
    icer.mn=mean(icer), 
    icer.md=median(icer),
    icer.lo95=quantile(icer,c(0.025)), 
    icer.hi95=quantile(icer,c(0.975)),
    vac_tp.mn=mean(vac_tp), 
    vac_tp.md=median(vac_tp),
    vac_tp.lo95=quantile(vac_tp,c(0.025)), 
    vac_tp.hi95=quantile(vac_tp,c(0.975))
  ),
  by=.(iso3, region, scen_id, with_SB, cet_name)
]
icers.countries <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][icers.countries, on=.(scen_id)] # add scenario names

####################      *** Prob. Cost-Effective ***        ##################

print("Probability Cost-Effective")

# global
ce.global <- icers.dt[
  ,.(nmb = sum(nmb)), 
  by=.(scen_id,with_SB,samp_id,cet_name)
][
  ,.(prob.CE = sum(nmb>0)/.N),
  by=.(scen_id,with_SB,cet_name)
] 
ce.global <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][ce.global, on=.(scen_id)] # add scenario names

# by region
ce.regions <- icers.dt[
  ,.(nmb = sum(nmb)), 
  by=.(region,scen_id,with_SB,cet_name,samp_id)
][
  ,.(prob.CE = sum(nmb>0)/.N),
  by=.(region,scen_id,with_SB,cet_name)
]
ce.regions <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][ce.regions, on=.(scen_id)] # add scenario names

# by country

ce.countries <- icers.dt[
  ,.(prob.CE = sum(nmb>0)/.N),
  by=.(iso3,region,scen_id,with_SB,cet_name)
]
ce.countries <- scen_pars[,.(scen_id,base_name,vac_name,scen_name)][ce.countries, on=.(scen_id)] # add scenario names

######################        *** WRITE RESULTS  ***        ####################

print("write output")

saveRDS(base.global,paste0(path.io,"base.global_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(base.regions,paste0(path.io,"base.regions_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(base.countries,paste0(path.io,"base.countries_",Sys.Date(),"_",reg_type,".rds"))

saveRDS(scenarios.global,paste0(path.io,"scenarios.global_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(scenarios.regions,paste0(path.io,"scenarios.regions_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(scenarios.countries,paste0(path.io,"scenarios.countries_",Sys.Date(),"_",reg_type,".rds"))

saveRDS(icers.global,paste0(path.io,"icers.global_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(icers.regions,paste0(path.io,"icers.regions_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(icers.countries,paste0(path.io,"icers.countries_",Sys.Date(),"_",reg_type,".rds"))

saveRDS(ce.global,paste0(path.io,"prob.CE.global_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(ce.regions,paste0(path.io,"prob.CE.regions_",Sys.Date(),"_",reg_type,".rds"))
saveRDS(ce.countries,paste0(path.io,"prob.CE.countries_",Sys.Date(),"_",reg_type,".rds"))

# report run-time
print(paste0("Done: ",Sys.time() - start_time))

######################     ***   TABLES & PLOTS    ***     #####################

# labels

if (reg_type=="sdg"){
  reg.lbl  <- c(
    'Central_and_Southern_Asia'="Central & Southern Asia",
    'Sub-Saharan_Africa'="Sub-Saharan Africa",
    'Europe_and_Northern_America'="Europe & Northern America",
    'Northern_Africa_Western_Asia'="Northern Africa & Western Asia",
    'Latin_America_and_Caribbean'="Latin America & Caribbean",
    'Oceania'="Oceania",
    'Eastern_and_South-Eastern_Asia'="Eastern & South Eastern Asia"
  ) } else if (reg_type=="wb") {
    reg.lbl  <- c(
      'Low_income'="Low income",
      'Lower_middle_income'="Lower middle income",
      'Upper_middle_income'="Upper middle income",
      'High_income'="High income"
    ) } else {print("Error: reg_type undefined")}

norm.lbl <- c(
  'least favourable'="Normative assumptions:\nleast favourable",
  'most favourable'="most favourable"
)
disc.lbl <- c('0'="0% QALY disc.",'0.03'="3% QALY disc.")
SB.lbl   <- c('No'="No SB QALYs",'Yes'="With SB QALYs")

# function to add label for conservative/optimistic scenarios

con_opt <- function(dt,s=scen_pars){
  dt <- merge(dt,s[,.(scen_id,disc.benefits)],by="scen_id")
  dt <- dt[
    disc.benefits == 0.03 & with_SB=="No" & cet_name=="Empirical",
    norm_scen := "least favourable"
  ][
    disc.benefits == 0 & with_SB=="Yes" & cet_name=="1xGDPPC",
    norm_scen := "most favourable"
  ]
  return(dt)
}

icers.global <- con_opt(icers.global,scen_pars)
icers.regions <- con_opt(icers.regions,scen_pars)
icers.countries <- con_opt(icers.countries,scen_pars)

ce.global <- con_opt(ce.global,scen_pars)
ce.regions <- con_opt(ce.regions,scen_pars)
ce.countries <- con_opt(ce.countries,scen_pars)

# add CET threshold
icers.countries <- dl$countries.dt[,.(iso3, gdppc = wb_gdp_per_capita_usd_2020)][icers.countries,on=.(iso3)]
scenarios.countries <- dl$countries.dt[,.(iso3, gdppc = wb_gdp_per_capita_usd_2020)][scenarios.countries,on=.(iso3)]
ce.countries <- dl$countries.dt[,.(iso3, gdppc = wb_gdp_per_capita_usd_2020)][ce.countries,on=.(iso3)]

######################         ***   TABLES   ***          #####################

my_format <- function(x) {prettyNum(signif(x,3),big.mark=",")}
# my_format <- function(x) {label_number_si(accuracy=3)(x)}
my_cis <- function(x,y,z,scale=1){
  sprintf("%s \n(%s, %s)",
          my_format(x/scale),
          my_format(y/scale),
          my_format(z/scale)
  )}

# add scenario display names

# add region names
# icers.countries <- dl$regions.dt[,.(iso3,region=sdg)][icers.countries, on=.(iso3)] 
# scenarios.countries <- dl$regions.dt[,.(iso3,region=sdg)][scenarios.countries, on=.(iso3)]
# ce.countries <- dl$regions.dt[,.(iso3,region=sdg)][ce.countries, on=.(iso3)]

# Base case vaccine impact + vaccine preterm impact scenario

# regions
tab.regions.impact <- rbindlist(list(
  dcast(scenarios.regions[scen_id==9 & name %in% c("doses")       
                   ,.(region = region,
                     Description = "Number of women vaccinated \n(millions)",
                     value = my_format(inc.md/10^6)
                   )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("cost.vac_prog")
                    ,.(Description = "Vaccine programme costs \n(discounted; $ millions)",
                       region = region,
                       value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                    )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("cost.acute_hospital")
                        ,.(Description = "Acute healthcare costs \n(discounted; $ millions)",
                           region = region,
                           value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                        )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("cost.longterm")
                        ,.(Description = "Long-term healthcare costs \n(discounted; $ millions)",
                           region = region,
                           value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                        )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("cost.total")
                        ,.(Description = "Total incremental costs \n(discounted; $ millions)",
                           region = region,
                           value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                        )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("eod")
                          ,.(Description = "EOGBS cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("lod")
                          ,.(Description = "LOGBS cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("modsev_ndi")
                          ,.(Description = "Moderate & severe NDI cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("gbs_deaths")
                          ,.(Description = "GBS deaths \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("gbs_stillbirth")
                          ,.(Description = "GBS stillbirths \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==11 & name %in% c("gbs_preterm")
                          ,.(Description = "GBS associated preterm births* \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(discounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==9 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(discounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==11 & name %in% c("qaly_loss.gbs_preterm")
                          ,.(Description = "QALYs from averted preterm births* \n(discounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==10 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(undiscounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==10 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(undiscounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==12 & name %in% c("qaly_loss.gbs_preterm")
                          ,.(Description = "QALYs from averted preterm births* \n(undiscounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  )
))

# global
tab.global.impact <- rbindlist(list(
  dcast(scenarios.global[scen_id==9 & name %in% c("doses")       
                          ,.(region = "Global",
                             Description = "Number of women vaccinated \n(millions)",
                             value = my_format(inc.md/10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("cost.vac_prog")
                          ,.(Description = "Vaccine programme costs \n(discounted; $ millions)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("cost.acute_hospital")
                          ,.(Description = "Acute healthcare costs \n(discounted; $ millions)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("cost.longterm")
                          ,.(Description = "Long-term healthcare costs \n(discounted; $ millions)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("cost.total")
                          ,.(Description = "Total incremental costs \n(discounted; $ millions)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("eod")
                          ,.(Description = "EOGBS cases \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("lod")
                          ,.(Description = "LOGBS cases \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("modsev_ndi")
                          ,.(Description = "Moderate & severe NDI cases \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("gbs_deaths")
                          ,.(Description = "GBS deaths \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("gbs_stillbirth")
                          ,.(Description = "GBS stillbirths \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==11 & name %in% c("gbs_preterm")
                          ,.(Description = "GBS associated preterm births* \n(thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(discounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==9 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(discounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==11 & name %in% c("qaly_loss.gbs_preterm")
                          ,.(Description = "QALYs from averted preterm births* \n(discounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==10 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(undiscounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==10 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(undiscounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==12 & name %in% c("qaly_loss.gbs_preterm")
                          ,.(Description = "QALYs from averted preterm births* \n(undiscounted; thousands)",
                             region = "Global",
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  )
))

tab.impact <- cbind(tab.regions.impact,tab.global.impact[,.("Global^" = Global)])
# tab.impact <- merge(tab.global.impact,tab.regions.impact,by="description")

# strip underscores from region names
names(tab.impact) <- gsub("_"," ",names(tab.impact))

write.csv(tab.impact,file=paste0(path.io,"Table_1-incremental_outcomes",Sys.Date(),"_",reg_type,".csv"))


# High impact vaccine scenario

# regions
tab.regions.impact <- rbindlist(list(
  dcast(scenarios.regions[scen_id==19 & name %in% c("doses")       
                          ,.(region = region,
                             Description = "Number of women vaccinated \n(millions)",
                             value = my_format(inc.md/10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("cost.vac_prog")
                          ,.(Description = "Vaccine programme costs \n(discounted; $ millions)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("cost.acute_hospital")
                          ,.(Description = "Acute healthcare costs \n(discounted; $ millions)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("cost.longterm")
                          ,.(Description = "Long-term healthcare costs \n(discounted; $ millions)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("cost.total")
                          ,.(Description = "Total incremental costs \n(discounted; $ millions)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("eod")
                          ,.(Description = "EOGBS cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("lod")
                          ,.(Description = "LOGBS cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("modsev_ndi")
                          ,.(Description = "Moderate & severe NDI cases \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("gbs_deaths")
                          ,.(Description = "GBS deaths \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("gbs_stillbirth")
                          ,.(Description = "GBS stillbirths \n(thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(discounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==19 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(discounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==20 & name %in% c("qaly_loss.gbs_disease")
                          ,.(Description = "QALYs from averted GBS disease \n(undiscounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.regions[scen_id==20 & name %in% c("qaly_loss.gbs_stillbirth")
                          ,.(Description = "QALYs from averted stillbirths \n(undiscounted; thousands)",
                             region = region,
                             value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                          )], Description ~ region, value.var = "value"
  )
))

# global
tab.global.impact <- rbindlist(list(
  dcast(scenarios.global[scen_id==19 & name %in% c("doses")       
                         ,.(region = "Global",
                            Description = "Number of women vaccinated \n(millions)",
                            value = my_format(inc.md/10^6)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("cost.vac_prog")
                         ,.(Description = "Vaccine programme costs \n(discounted; $ millions)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("cost.acute_hospital")
                         ,.(Description = "Acute healthcare costs \n(discounted; $ millions)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("cost.longterm")
                         ,.(Description = "Long-term healthcare costs \n(discounted; $ millions)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("cost.total")
                         ,.(Description = "Total incremental costs \n(discounted; $ millions)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^6)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("eod")
                         ,.(Description = "EOGBS cases \n(thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("lod")
                         ,.(Description = "LOGBS cases \n(thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("modsev_ndi")
                         ,.(Description = "Moderate & severe NDI cases \n(thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("gbs_deaths")
                         ,.(Description = "GBS deaths \n(thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("gbs_stillbirth")
                         ,.(Description = "GBS stillbirths \n(thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.lo95,inc.hi95,scale=10^3)
                         )], Description ~ region, value.var = "value"
  ),
 
  dcast(scenarios.global[scen_id==19 & name %in% c("qaly_loss.gbs_disease")
                         ,.(Description = "QALYs from averted GBS disease \n(discounted; thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==19 & name %in% c("qaly_loss.gbs_stillbirth")
                         ,.(Description = "QALYs from averted stillbirths \n(discounted; thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==20 & name %in% c("qaly_loss.gbs_disease")
                         ,.(Description = "QALYs from averted GBS disease \n(undiscounted; thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                         )], Description ~ region, value.var = "value"
  ),
  dcast(scenarios.global[scen_id==20 & name %in% c("qaly_loss.gbs_stillbirth")
                         ,.(Description = "QALYs from averted stillbirths \n(undiscounted; thousands)",
                            region = "Global",
                            value = my_cis(inc.md,inc.hi95,inc.lo95,scale=-10^3)
                         )], Description ~ region, value.var = "value"
  )
))

tab.impact <- cbind(tab.regions.impact,tab.global.impact[,.("Global^" = Global)])

# strip underscores from region names
names(tab.impact) <- gsub("_"," ",names(tab.impact))

write.csv(tab.impact,file=paste0(path.io,"Supplmentary_table_1-high-coverage_incremental_outcomes",Sys.Date(),"_",reg_type,".csv"))


######################         ***   FIGURES  ***          #####################

# *** Figures for runs with SDG regions ***

if (reg_type=="sdg") { 

# FIGURE 2a - Base case global NMB under different normative assumptions

fig2a <- ggplot(
  data = icers.global[scen_name=="Base case"]
) +
  geom_pointrange(
    aes(
      x = cet_name, 
      y=nmb.md, 
      ymin=nmb.lo95, 
      ymax=nmb.hi95, 
      color = cet_name,
      # group = interaction(cet_name,with_SB,disc.benefits)
      ),
    position = position_dodge(width=0.5)
  ) +
  scale_color_brewer(palette="Set1") +
  facet_grid(.~disc.benefits+with_SB,
             labeller = labeller(
               disc.benefits = disc.lbl,
               with_SB = SB.lbl
             )) +
  scale_y_continuous(labels = scales::label_number_si()) +
  labs(
    # title="Base case global Net Monetary Benefit of GBS vaccination\n under different normative assumptions",
    title="A",
    color="CET",
    # caption="GDPPC = GDP per capita; SB = stillbirth; \nCET = Cost-effectiveness Threshold"
  ) +
  ylab("Net Monetary Benefit (USD)") +
  xlab("") +
  theme_minimal() +
  theme(legend.position="right",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  NULL

# ggsave(
#   paste0(path.io,"FIG2a-nmb_global_",Sys.Date(),"_",reg_type,".png"), 
#   fig2a, width = 7, height = 5, units = "in"
# )

# FIGURE 2b - Regional NMB under base case assumptions

fig2b <- ggplot(
  data = icers.regions[scen_name=="Base case" & !is.na(norm_scen)]
) +
  geom_pointrange(
    aes(
      x = region, 
      y=nmb.md, 
      ymin=nmb.lo95, 
      ymax=nmb.hi95, 
      color = region,
    ),
    position = position_dodge(width=0.5)
  ) +
  facet_wrap(
    .~norm_scen,
    labeller = labeller(norm_scen = norm.lbl),
    ncol=2#,
    # scales="free"
  ) +
  scale_y_continuous(
    trans=ggallin::ssqrt_trans,breaks=(c(-0.5,0,0.5,2,5,10)*10^9),
    labels = scales::label_number_si()
  ) +
  theme(legend.position="top") +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  labs(
    # title="Base case regional Net Monetary Benefit of GBS vaccination",
    title="B",
    # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  ylab("Net Monetary Benefit (USD)") +
  xlab("") +
  theme_minimal() +
  theme(legend.position="right",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  guides(color=guide_legend(ncol=1, bycol=TRUE)) +
  NULL

# ggsave(
#   paste0(path.io,"FIG2b-nmb_regions_base_",Sys.Date(),"_",reg_type,".png"), 
#   fig2b, width = 7, height = 5, units = "in"
# )

# FIGURE 2 - Combine Fig 2a and b

fig2 <- grid.arrange(fig2a, fig2b, ncol=1, nrow=2, heights=c(1,1))
ggsave(  
  paste0(path.io,"FIG2-nmb_base_combined_",Sys.Date(),"_",reg_type,".png"),
  fig2, width=7, height=10, units="in"
)

# FIGURE 3 - base case probability cost effective by country

fig3 <- ggplot(
  data=ce.countries[scen_name=="Base case" & !is.na(norm_scen)],
) +
  geom_point(aes(x=gdppc, y=prob.CE, color=region), size=1, alpha=0.6) +
  facet_wrap(
    .~norm_scen,
    ncol = 2,
    labeller = labeller(norm_scen = norm.lbl) 
  ) +
  scale_x_continuous(trans="log10") +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  labs(
    # title="Base case probability vaccine is cost-effective by country",
    # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  ylab("Probability cost-effective") +
  xlab("Countries ordered by GDP per capita") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  guides(color=guide_legend(ncol=3, bycol=TRUE)) +
  NULL

ggsave(
  paste0(path.io,"FIG3-p_CE_base_",Sys.Date(),"_",reg_type,".png"), 
  fig3, width = 7, height = 5, units = "in"
)

# FIGURE 4 - Global NMB scenario analysis

fig4 <- ggplot(
  data = icers.global[!is.na(norm_scen)]
) +
  geom_pointrange(
    aes(
      y = reorder(scen_name,scen_id), 
      x=nmb.md, 
      xmin=nmb.lo95, 
      xmax=nmb.hi95,
      # shape = norm_scen,
      color = norm_scen
    ),
    position = position_dodge(width=0.5)
  ) +
  scale_color_brewer(palette="Set1",
    name="Normative assumptions"
  ) +
  scale_x_continuous(
    trans=ggallin::ssqrt_trans,breaks=(c(-1.0000001,0,1.0000001,5,10,20,30,40,50)*1e9),
    labels = scales::label_number_si()
    ) +
  labs(
    # title="Regional Net Monetary Benefit of GBS vaccination \nunder different scenarios",
    shape="normative assumptions",
    # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  xlab("Net Monetary Benefit (USD)") +
  ylab("Scenario") +
  theme_minimal() +
  theme(legend.position="bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        # strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  NULL

ggsave(
  paste0(path.io,"FIG4-nmb_global_scen_",Sys.Date(),"_",reg_type,".png"), 
  fig4, width = 7, height = 5, units = "in"
)

# SUPPLEMENTARY FIGURE 1 - base case probability cost effective under different normative assumptions

sfig1 <- ggplot(
  data=icers.regions[scen_name=="Base case"],
) +
  geom_pointrange(
    aes(
      x = region, 
      y=nmb.md, 
      ymin=nmb.lo95, 
      ymax=nmb.hi95, 
      color = region,
    ))+
  facet_grid(disc.benefits~cet_name+with_SB,
             labeller = labeller(
               disc.benefits = disc.lbl,
               with_SB = SB.lbl
             )) +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  scale_y_continuous(
    trans=ggallin::ssqrt_trans,breaks=(c(-0.5,0,0.5,2,5,10)*10^9),
    labels = scales::label_number_si()
  ) +
  labs(
    # title="Regional NMB for vaccine base case \nunder different normative assumptions",
    # caption="GDPPC = GDP per capita; SB = stillbirth"
  ) +
  ylab("Net Monetary Benefit (USD)") +
  xlab("") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  guides(color=guide_legend(ncol=3, bycol=TRUE)) +
  NULL

ggsave(
  paste0(path.io,"SUP_FIG1-nmb_regions_norm_scen_",Sys.Date(),"_",reg_type,".png"), 
  sfig1, width = 7, height = 10, units = "in"
)

# SUPPLEMENTARY FIGURE 2 - Regional NMB scenario analysis

sfig2 <- ggplot(
  data = icers.regions[!is.na(norm_scen)]
) +
  geom_pointrange(
    aes(
      x = region, 
      y=nmb.md, 
      ymin=nmb.lo95, 
      ymax=nmb.hi95, 
      color = region,
      shape = norm_scen
    ),
    position = position_dodge(width=0.5)
  ) +
  facet_wrap(.~reorder(scen_name,scen_id),
             ncol=5,
             labeller = label_wrap_gen(width = 12)
  ) +
  scale_y_continuous(trans=ggallin::ssqrt_trans,breaks=(c(-0.5,0,0.5,2,5,10,20)*10^9),
    labels = scales::label_number_si()) +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  labs(
    # title="Regional Net Monetary Benefit under different GBS vaccination scenarios",
    shape="normative assumptions",
    # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  ylab("Net Monetary Benefit (USD)") +
  xlab("") +
  theme_minimal() +
  theme(legend.position="bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  NULL

ggsave(
  paste0(path.io,"SUP_FIG2-nmb_regions_scen_",Sys.Date(),"_",reg_type,".png"), 
  sfig2, width = 7, height = 10, units = "in"
)

# SUPPLEMENTARY FIGURE 3 - base case probability cost effective under different normative assumptions

sfig3 <- ggplot(
  data=ce.countries[scen_name=="Base case"],
) +
  geom_point(aes(x=gdppc, y=prob.CE, color=region), size=1, alpha=0.6) +
  facet_grid(disc.benefits~cet_name+with_SB,
             labeller = labeller(
               disc.benefits = disc.lbl,
               with_SB = SB.lbl
             )) +
  scale_x_continuous(trans="log10") +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  labs(
    # title="Probability vaccine basecase is cost-effective by country \nunder different normative assumptions",
    # caption="GDPPC = GDP per capita; SB = stillbirth"
  ) +
  ylab("Probability cost-effective") +
  xlab("Countries ordered by GDP per capita") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  guides(color=guide_legend(ncol=3, bycol=TRUE)) +
  NULL

ggsave(
  paste0(path.io,"SUP_FIG3-p_CE_norm_scen_",Sys.Date(),"_",reg_type,".png"), 
  sfig3, width = 7, height = 6, units = "in"
)

# SUPPLEMENTARY FIGURE 4 - base case probability cost effective under different vaccine scenarios

sfig4 <- ggplot(
  data=ce.countries[!is.na(norm_scen)],
) +
  geom_point(aes(x=gdppc, y=prob.CE, color=region), size=1, alpha=0.6) +
  facet_grid(
    reorder(scen_name,scen_id)~norm_scen,
    labeller = labeller(norm_scen = norm.lbl, .rows = label_wrap_gen(width = 11))
  ) +
  scale_x_continuous(trans="log10") +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  labs(
    # title = "Probability vaccine is cost-effective by country under different scenarios",
    # caption = "least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  ylab("Probability cost-effective") +
  xlab("Countries ordered by GDP per capita") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  guides(color=guide_legend(ncol=3, bycol=TRUE)) +
  NULL

ggsave(
  paste0(path.io,"SUP_FIG4-p_CE_scen_",Sys.Date(),"_",reg_type,".png"), 
  sfig4, width = 7, height = 10, units = "in"
)

# SUPPLEMENTARY FIGURE 5 - SDG Regional base case vaccine threshold prices

sfig5 <- ggplot(
  data=icers.countries[scen_name=="Base case" & !is.na(norm_scen)]
) +
  geom_boxplot(aes(
    x=region, 
    y=vac_tp.md, 
    color=region),
    outlier.shape = NA
  ) +
  facet_wrap(
    .~norm_scen,
    ncol = 2,
    # scales = "free",
    labeller = labeller(norm_scen = norm.lbl) 
  ) +
  scale_color_brewer(palette="Set1",
    name="SDG Region",
    labels = reg.lbl
  ) +
  scale_y_continuous(trans=ggallin::ssqrt_trans,breaks=c(-10,-1,0,1,10,50,100,200,400,800)) +
  labs(
    # title="Distribution of threshold prices for a GBS vaccine by region",
    # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
  ) +
  ylab("Threshold vaccine price (USD)") +
  xlab("") +
  theme_minimal() +
  theme(legend.position="bottom",
        plot.background=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_rect(colour= "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour = "gray", size=0.5),
        panel.grid.minor = element_line(colour = "gray", size=0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(hjust=0, vjust=0),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust=0,size=9)
  ) +
  coord_cartesian(ylim=c(-10,800)) +
  guides(color=guide_legend(ncol=3, bycol=TRUE)) +
  NULL

ggsave(
  paste0(path.io,"SUP_FIG5-vac_tp_base_",Sys.Date(),"_",reg_type,".png"), 
  sfig5, width = 7, height = 6, units = "in"
)

} # end of figures for sdg region runs

# *** Figures for WB region runs ***

if (reg_type=="wb") {
  
  # reorder regions
  icers.regions$region <- factor(icers.regions$region,levels=c("Low_income","Lower_middle_income","Upper_middle_income","High_income"))
  icers.countries$region <- factor(icers.countries$region,levels=c("Low_income","Lower_middle_income","Upper_middle_income","High_income"))
  
  # FIGURE 5 - WB Regional base case vaccine threshold prices
  
  fig5 <- ggplot(
    data=icers.countries[scen_name=="Base case" & !is.na(norm_scen)]
  ) +
    geom_boxplot(aes(
      x=region, 
      y=vac_tp.md, 
      color=region),
      outlier.shape = NA
    ) +
    facet_wrap(
      .~norm_scen,
      ncol = 2,
      # scales = "free",
      labeller = labeller(norm_scen = norm.lbl) 
    ) +
    scale_color_brewer(palette="Set1",
      name="WB Income Group",
      labels = reg.lbl
    ) +
    scale_y_continuous(trans=ggallin::ssqrt_trans,breaks=c(-10,-1,0,1,10,50,100,200,400,800)) +
    labs(
      # title="Distribution of threshold prices for a GBS vaccine by World Bank income group",
      # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
    ) +
    ylab("Threshold vaccine price (USD)") +
    xlab("") +
    theme_minimal() +
    theme(legend.position="bottom",
          plot.background=element_rect(colour = "black", fill=NA, size=1),
          panel.border=element_rect(colour= "black", fill=NA, size=0.5),
          panel.grid.major = element_line(colour = "gray", size=0.5),
          # panel.grid.minor = element_line(colour = "gray", size=0.5),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text.x = element_text(hjust=0, vjust=0),
          plot.caption.position = "plot",
          plot.caption = element_text(hjust=0,size=9)
    ) +
    coord_cartesian(ylim=c(-10,800)) +
    NULL
  
  ggsave(
    paste0(path.io,"FIG5-vac_tp_base_",Sys.Date(),"_",reg_type,".png"), 
    fig5, width = 7, height = 6, units = "in"
  )
  
  # SUPPLEMENTARY FIGURE 6 - WB Regional base case vaccine threshold price scenario analysis
  
 sfig6 <- ggplot(
    data=icers.countries[!is.na(norm_scen)]
  ) +
    geom_boxplot(aes(
      y=region, 
      x=vac_tp.md, 
      color=region),
      outlier.shape = NA
    ) +
    facet_grid(
      reorder(scen_name, scen_id)~norm_scen,
      labeller = labeller(norm_scen = norm.lbl, .rows = label_wrap_gen(width = 11))
    ) +
   scale_color_brewer(palette="Set1",
      name="WB Income Group",
      labels = reg.lbl
    ) +
   scale_x_continuous(trans=ggallin::ssqrt_trans,breaks=c(-10,-1,0,1,10,50,100,200,400,800)) +
    labs(
      # title="Distribution of regional vaccine threshold prices by World Bank income group \nunder different vaccine scenarios",
      # caption="least favourable assumptions: empirical CET, 3% QALY discounting, excludes stillbirth QALYs\nmost favourable assumptions: 1 x GDP per capita CET, 0% QALY discounting, includes stillbirth QALYs"
    ) +
    xlab("Threshold vaccine price (USD)") +
    ylab("") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.background=element_rect(colour = "black", fill=NA, size=1),
          panel.border=element_rect(colour = "black", fill=NA, size=0.5),
          panel.grid.major = element_line(colour = "gray", size=0.5),
          panel.grid.minor = element_line(colour = "gray", size=0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x = element_text(hjust=0, vjust=0),
          plot.caption.position = "plot",
          plot.caption = element_text(hjust=0,size=9)
    ) +
    guides(color=guide_legend(ncol=2, bycol=TRUE)) +
    coord_cartesian(xlim=c(-50,800)) +
    NULL
  
  ggsave(
    paste0(path.io,"SUP_FIG6-vac_tp_scen_",Sys.Date(),"_",reg_type,".png"), 
    sfig6, width =7, height = 10, units = "in"
  )
  
}

} # END OF RUN FOR BOTH REGION GROUPINGS








