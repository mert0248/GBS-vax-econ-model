################ *** FUNCTION TO BUILD & EVALUATE TREE  ***    #################

require(data.table)

eval_tree <- function(pars){
  
  # parameters as a data table
  dt <- as.data.table(pars)
  
  # *** Risks ***
  # colonisation
  dt[, colonised      := p.mat_col]
  
  # select vaccine coverage assumption
  if (pars$vac_cov_assumption==0){
    dt[, p.vaccine_coverage := 0]
  } else if (pars$vac_cov_assumption==1){
    dt[, p.vaccine_coverage := p.anc1]
  } else { # default to ANC4
    dt[, p.vaccine_coverage := p.anc4]
  }
  
  if (pars$n.doses==2) { # 2-dose regimen
  
    # vaccine coverage is adjusted by p.dose2 -- the probability that a mother receiving the first dose also receives the second
    # assumes no efficacy from a single dose

    # stillbirth
    dt[, gbs_stillbirth := colonised * (risk.stillbirth.col - risk.stillbirth.no_col) * (1 - p.ve_stillbirth * p.vaccine_coverage * p.dose2)]
    dt[, all_stillbirth := risk.stillbirth.no_col + gbs_stillbirth]
    dt[, all_livebirth  := 1 - all_stillbirth]
    
    # preterm births - vaccine effective coverage is adjusted for p.preterm_covered the estimated proportion of preterm births that occur after vaccine is given
    dt[, gbs_preterm    := (colonised * (1 - risk.stillbirth.no_col) - gbs_stillbirth) * (risk.preterm.col - risk.preterm.no_col) * (1 - p.ve_preterm * p.vaccine_coverage * p.preterm_covered * p.dose2)]
    dt[, all_preterm    := all_livebirth * risk.preterm.no_col + gbs_preterm]
    dt[, all_term       := all_livebirth - all_preterm]
    
    # early-onset disease: live births to colonised mothers x risk if colonised x vaccine effect (NB including extra live births from averted gbs stillbirths)
    dt[, eod            := (colonised * (1 - risk.stillbirth.no_col) - gbs_stillbirth) * p.eod.col * (1 - p.ve_gbs_disease * p.vaccine_coverage * p.dose2)] 
  
  } else { # 1-dose regimen and non-vaccine scenarios
    
    # stillbirth
    dt[, gbs_stillbirth := colonised * (risk.stillbirth.col - risk.stillbirth.no_col) * (1 - p.ve_stillbirth * p.vaccine_coverage)]
    dt[, all_stillbirth := risk.stillbirth.no_col + gbs_stillbirth]
    dt[, all_livebirth  := 1 - all_stillbirth]
    
    # preterm births - vaccine effective coverage is adjusted for p.preterm_covered the estimated proportion of preterm births that occur after vaccine is given
    dt[, gbs_preterm    := (colonised * (1 - risk.stillbirth.no_col) - gbs_stillbirth) * (risk.preterm.col - risk.preterm.no_col) * (1 - p.ve_preterm * p.vaccine_coverage * p.preterm_covered)]
    dt[, all_preterm    := all_livebirth * risk.preterm.no_col + gbs_preterm]
    dt[, all_term       := all_livebirth - all_preterm]
    
    # early-onset disease: live births to colonised mothers x risk if colonised x vaccine effect (NB including extra live births from averted gbs stillbirths)
    dt[, eod            := (colonised * (1 - risk.stillbirth.no_col) - gbs_stillbirth) * p.eod.col * (1 - p.ve_gbs_disease * p.vaccine_coverage)] 
    
  }
  
  # late-onset disease: note using ratio so any vax effect applied to eod is already applied here
  dt[, lod            := eod * prop.lod / (1 - prop.lod)]
  # deaths due to gbs disease - adjusted for coverage of skilled birth attendance
  dt[, eod_deaths      := eod * (coverage_SBA * p.gbs_death.eod + (1 - coverage_SBA) * p.eod_death_no_sba)]
  dt[, lod_deaths      := lod * p.gbs_death.lod]
  dt[, gbs_deaths      := eod_deaths + lod_deaths]
  # meningitis and sepsis
  dt[, meningitis     := (eod - eod_deaths) * p.men.eod + (lod - lod_deaths) * p.men.lod]
  dt[, sepsis         := (eod - eod_deaths) * (1 - p.men.eod) + (lod - lod_deaths) * (1 - p.men.lod)]
  
  # excess neuro-developmental impairment
  dt[, any_ndi        := meningitis * (1 - p.no_ndi.men) + sepsis * (1 - p.no_ndi.sep)]
  dt[, mild_ndi       := meningitis * p.mild_ndi.men + sepsis * p.mild_ndi.sep]
  dt[, mod_ndi        := meningitis * p.mod_ndi.men  + sepsis * p.mod_ndi.sep]
  dt[, sev_ndi        := meningitis * p.sev_ndi.men  + sepsis * p.sev_ndi.sep]
  dt[, modsev_ndi     := mod_ndi + sev_ndi]
  # absolute neuro-developmental impairment
  dt[, abs_modsev_ndi := meningitis * p.modsev_ndi_abs.men  + sepsis * p.modsev_ndi_abs.sep]
  dt[, abs_any_ndi    := meningitis * p.any_ndi_abs.men  + sepsis * p.any_ndi_abs.sep]
  
  # *** QALYS *** *** 
  
  # life expectancy
  le <- floor(pars[["life_expectancy"]])
  # discounted life years
  rate <- pars[["disc.benefits"]]
  ly <- sum((1+rate)^-seq(1:le))
 
  # QALY losses for GBS associated outcomes
  
  # qaly loss for child acute illness
  dt[, qaly_loss.sep      := sepsis     * u.sep * acute_los/365] # decrements for hospitalisation due to acute sepsis
  dt[, qaly_loss.men      := meningitis * u.men * acute_los/365] # decrements for hospitalisation due to acute meningitis
  
  # qaly loss for deaths - assumes deaths occur proportionately across term/preterm branches
  dt[, qaly_loss.deaths   := gbs_deaths * ly * ((1-u.term) * all_preterm/all_livebirth + (1-u.preterm) * all_term/all_livebirth)]
  
  # qaly loss for NDI - *** currently assumed additive with preterm qaly loss and no change in life expectancy in survivors ***
  dt[, qaly_loss.mild_ndi := mild_ndi * ly * u.mild_ndi]
  dt[, qaly_loss.mod_ndi  := mod_ndi  * ly * u.mod_ndi]
  dt[, qaly_loss.sev_ndi  := sev_ndi  * ly * u.sev_ndi]
  dt[, qaly_loss.any_ndi  := qaly_loss.mild_ndi + qaly_loss.mod_ndi + qaly_loss.sev_ndi]
  
  dt[, qaly_loss.gbs_disease  := qaly_loss.sep + qaly_loss.men + qaly_loss.deaths + qaly_loss.mild_ndi + qaly_loss.mod_ndi + qaly_loss.sev_ndi]
  
  dt[, qaly_loss.gbs_preterm    := gbs_preterm * ly * u.preterm]
  # assume the fraction of averted gbs_stillbirths born preterm are the same as the general population
  dt[, qaly_loss.gbs_stillbirth := gbs_stillbirth * ly * ((1-u.term) * all_preterm/all_livebirth + (1-u.preterm) * all_term/all_livebirth)]
  dt[, qaly_loss.gbs_with_stillbirth  := qaly_loss.gbs_disease + qaly_loss.gbs_stillbirth]
  
  # qalys for live births before losses due to gbs 
  # dt[, qalys.term     := all_term    * (1 - u.term)    * ly] # qalys for all live term births 
  # dt[, qalys.preterm  := all_preterm * (1 - u.preterm) * ly] # qalys for all live preterm births
  
  # total qaly loss
  dt[, qaly_loss.total    := qaly_loss.gbs_with_stillbirth + qaly_loss.gbs_preterm]

  # *** Costs ***
  # vaccine program costs
  
  if (pars$n.doses==2) { # 2-dose regimen
    dt[, doses:= p.vaccine_coverage * (1 + p.dose2)]
  } else {
    dt[, doses:= p.vaccine_coverage]
  }
  dt[, cost.vac_delivery := doses * c.vac_delivery * vac_delivery_cost_mult] # vac_delivery_cost_mult used to scale delivery costs in scenario analysis
  dt[, cost.vac_price    := doses * c.vac_price]
  dt[, cost.vac_prog     := cost.vac_delivery + cost.vac_price]
  
  # acute hospitalisation costs
  dt[, cost.acute_hospital := (eod + lod) * c.acute_hospital ]
  
  # discounted long-term costs over lifetime
  rate <- pars[["disc.costs"]]
  ly.c <- sum((1+rate)^-seq(1:le))  
  
  # applied to mod/sev ndi based on proportion of acute costs
  dt[, cost.longterm := (mod_ndi + sev_ndi) * ly.c * c.mod_sev_ndi_as_prop_acute * c.acute_hospital]
  dt[, cost.h_care.total := cost.acute_hospital + cost.longterm]
  
  # total cost
  dt[, cost.total := cost.vac_prog + cost.acute_hospital + cost.longterm]

       
  # *** Fully vaccinated***
  if (pars$n.doses==2) { # 2-dose regimen
    dt[, fully_vaccinated := p.vaccine_coverage * p.dose2]
  } else {
    dt[, fully_vaccinated := p.vaccine_coverage]
  }
  
  # add sample id
  dt[, samp_id := .I]
  
  # return selected outputs
  dt <- dt[,.(
    # colonised,
    gbs_stillbirth,
    # all_stillbirth,
    # all_livebirth,
    # all_term,
    gbs_preterm,
    # all_preterm,
    eod,
    lod,
    meningitis,
    sepsis,
    # eod_deaths,
    # lod_deaths,
    gbs_deaths,
    any_ndi,
    # mild_ndi,
    # mod_ndi,
    # sev_ndi,
    modsev_ndi,
    # abs_modsev_ndi,
    # abs_any_ndi,
    # qalys.term,
    # qalys.preterm,
    qaly_loss.gbs_disease,
    qaly_loss.gbs_preterm,
    qaly_loss.gbs_stillbirth,
    # qaly_loss.gbs_with_stillbirth, 
    # qaly_loss.sep,
    # qaly_loss.men,
    # qaly_loss.deaths,
    qaly_loss.any_ndi,
    # qaly_loss.mild_ndi,
    # qaly_loss.mod_ndi,
    # qaly_loss.sev_ndi,
    qaly_loss.total,
    cost.vac_delivery,
    cost.vac_price,
    cost.vac_prog,
    cost.acute_hospital,
    cost.longterm,
    cost.h_care.total,
    cost.total,
    doses,
    fully_vaccinated,
    samp_id
  )]
  
  return (dt)  
}