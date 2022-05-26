##########################        ***    HELPER FUNCTIONS   ***         ############################

###############################          FUNCTION: OR2R            #################################
# Calculate the corresponding risk in the exposed and unexposed groups, given
#
# i  - the incidence of an outcome 
# p  - the prevalence of a risk factor and
# or - the odds ratio of the outcome given the exposure

or2r <- function(i,p,or){
  
  
  # find quadratic root
  a = 1-or
  b = 1-p-i+or*(p+i)
  c = -or*p*i
  
  s1 = (-b + (b^2 - 4*a*c)^0.5)/(2*a)
  #s2 = (-b - (b^2 - 4*a*c)^0.5)/(2*a) # I think we can just ignore this root
  
  re = s1/p         # risk in exposed
  r0 = (i-s1)/(1-p) # risk in unexposed
  
  return(list(re=re,r0=r0))
}


############################         FUNCTION: GET PARAMETERS        ###############################
#reads in model parameter from data.parameters and for sthochastic runs draws random samples 
#using LHS otherwise returns a list of fixed values

get_samples <- function(data.parameters,stochastic,n_samp,n_dist){
  
  require("lhs") #support for latin hypercube sampling

  #generate random latin hypercube for n_samp samples and n_dist parameter distribution
  latin.sample <- randomLHS(n_samp,n_dist)
  
  #generate global set of parameter distributions
  param.names <- data.parameters[,full_name] #get list of parameters from data.parameters
  param.list <- vector("list",length(param.names)) #list to store parameter distributions
  names(param.list) = param.names
  
  for (i in 1:length(param.names)){
    param.name <- param.names[i]
    if (stochastic){ #if stochastic then draws samples from appropriate distributions
      param.list[[i]] <- with(data.parameters[full_name == param.name, ],
                              apply_latin_dist(latin.sample[,i],
                                               dist_name, dist.1, dist.2, value)) 
    } else { #otherwise return single fixed value
      param.list[[i]] = data.parameters[full_name == param.name, value]
    }
  }
  return(param.list) #return list of sampled parameter values
}

################################      FUNCTION: APPLY DIST        ##################################
# calls appropriate distribution function based on passed in parameters
# n is number of sample draws, dist is the distribution name,
# dist.1 and dist.2 are paramters used to define a distribution of type dist
# currently handles only fixed, uniform, normal, beta and gamma distributions

apply_dist <- function(n,dist,dist.1,dist.2,value){
  if (dist == "fixed"){
    r = value
  } else if (dist == "uniform"){
    r = runif(n,dist.1,dist.2)
  } else if(dist == "normal"){
    r = rnorm(n,dist.1,dist.2)
  } else if (dist == "beta"){
    r = rbeta(n,dist.1,dist.2)
  } else if(dist == "lognormal"){
    r = rlnorm(n,dist.1,dist.2)
  } else if(dist == "gamma"){
    r = rgamma(n,shape=dist.1,rate=dist.2) #using shape/rate
  } else {
    print("error in apply_dist() - distribution name not defined")
    r = NULL #if distribution name not supported return NULL
  }
  return(r)
}

##################################     FUNCTION: APPLY LATIN DIST      ######################################
# calls appropriate distribution function based on passed in parameters 
# k is vector a vector of samples to draw based on LHS, dist is the distribution name, 
# dist.1 and dist.2 are paramters used to define a distribution of type dist 
# currently handles only fixed, uniform, normal, beta and gamma distributions

apply_latin_dist <- function(k,dist,dist.1,dist.2,value){
  if (dist == "fixed"){
    r = value
  } else if  (dist == "uniform"){
    r = qunif(k,dist.1,dist.2)
  } else if(dist == "normal"){
    r = qnorm(k,dist.1,dist.2)
  } else if (dist == "beta"){
    r = qbeta(k,dist.1,dist.2)
  } else if (dist == "lognormal"){
    r = qlnorm(k,dist.1,dist.2) 
  } else if (dist == "gamma"){
    r = qgamma(k,shape=dist.1,rate=dist.2) # using shape/rate
  } else {
    print("error in apply_latin_dist() - distribution name not defined")
    r = NULL # if distribution name not supported return NULL
  }
  return(r)
}


