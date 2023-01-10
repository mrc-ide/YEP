# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YEP package
#-------------------------------------------------------------------------------
#' @title MCMC
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param log_params_ini Initial values of parameters to be varied These should always be the log() values of the
#'   actual parameters, ordered as follows:
#'   1) Parameters controlling the value of spillover force of infection FOI, either a) a number of FOI values equal
#'   to the total number of regions to be considered or b) a number of environmental coefficients used to calculate
#'   FOI values from environmental covariates equal to the number of environmental covariates listed in the
#'   enviro_data frame. Values should be in alphabetical order by region in case (a) or in the order of the columns
#'   in the environmental data frame in case (b).
#'   2) If the basic reproduction number for human-human transmission R0 is to be varied (i.e. type is set to
#'   "FOI+R0" or "FOI+R0 enviro"), parameters controlling the value of R0, either a) a number of R0 values equal to
#'   the total number of regions to be considered or b) a number of environmental coefficients used to calculate R0
#'   values from environmental covariates equal to the number of environmental covariates listed in the enviro_data
#'   frame. Values should be in alphabetical order by region in case (a) or in the order of the columns in the
#'   environmental data frame in case (b).
#'   3) Values of the additional parameters (vaccine efficacy vaccine_efficacy, severe case reporting probability
#'   p_rep_severe and fatal case reporting probability p_rep_death) if these are to be varied, in the order
#'   vaccine_efficacy->p_rep_severe->p_rep_death. If these parameters are to be varied, the values separately
#'   supplied to this function (see below) should be set to NULL, the default.
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param filename_prefix Prefix of names for output files
#' @param Niter Total number of steps to run
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param log_params_min Lower limits of varied parameter values if specified
#' @param log_params_max Upper limits of varied parameter values if specified
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param dt time increment in days (must be 1 or 5)
#' @param enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
#'   enviro"); set to NULL if not in use
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#' @param m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)
#' '
#' @export
#'
MCMC <- function(log_params_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                 filename_prefix="Chain",Niter=1,type=NULL,log_params_min=c(),log_params_max=c(),mode_start=0,
                 prior_type="zero",dt=1.0,enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,
                 p_rep_severe=NULL,p_rep_death=NULL,m_FOI_Brazil=1.0){

  #Check that initial, minimum and maximum parameters are in vectors of same sizes
  n_params=length(log_params_ini)
  assert_that(length(log_params_min)==n_params)
  assert_that(length(log_params_max)==n_params)

  extra_params=c()
  if(is.null(vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vaccine_efficacy")}
  if(is.null(p_rep_severe)==TRUE){extra_params=append(extra_params,"p_rep_severe")}
  if(is.null(p_rep_death)==TRUE){extra_params=append(extra_params,"p_rep_death")}
  if(is.null(m_FOI_Brazil)==TRUE){extra_params=append(extra_params,"m_FOI_Brazil")}

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data=input_data_process(input_data,obs_sero_data,obs_case_data)
  regions=names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions=length(regions)
  if(is.null(enviro_data)==FALSE){
    for(region in regions){assert_that(region %in% enviro_data$region)}
    enviro_data=subset(enviro_data,enviro_data$region %in% regions)
    }

  #Label parameters according to order and fitting type
  param_names=create_param_labels(type,input_data,enviro_data,extra_params)
  names(log_params_ini)=names(log_params_min)=names(log_params_max)=param_names

  #Run checks to ensure that number of parameters is correct for fitting type and number of regions/environmental
  #covariates
  checks<-mcmc_checks(log_params_ini,n_regions,type,log_params_min,log_params_max,prior_type,enviro_data,
                       R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death,m_FOI_Brazil)

  #Set up list of invariant parameter values to supply to other functions
  const_list=list(type=type,log_params_min=log_params_min,log_params_max=log_params_max,
                  mode_start=mode_start,prior_type=prior_type,dt=dt,enviro_data=enviro_data,
                  R0_fixed_values=R0_fixed_values,vaccine_efficacy=vaccine_efficacy,p_rep_severe=p_rep_severe,
                  p_rep_death=p_rep_death,m_FOI_Brazil=m_FOI_Brazil)

  ### find posterior probability at start ###
  out = MCMC_step(log_params=log_params_ini,input_data,obs_sero_data,obs_case_data,
                   chain_cov=1,adapt=0,like_current=-Inf,const_list)

  #MCMC setup
  chain=chain_prop=posterior_current=posterior_prop=flag_accept=chain_cov_all=NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  chain_cov=1

  #Iterative fitting
  for (iter in 1:Niter){
    #save current step
    log_params = out$log_params
    log_params_prop=out$log_params_prop
    like_current = out$like_current
    like_prop=out$like_prop
    accept = out$accept
    chain = rbind(chain, log_params)
    chain_prop=rbind(chain_prop,log_params_prop)
    posterior_current=rbind(posterior_current,like_current)
    posterior_prop=rbind(posterior_prop,like_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all,max(chain_cov))

    #Set output headings
    if(iter==1){
      colnames(chain)=colnames(chain_prop)=names(log_params_ini)
      for(i in 1:n_params){colnames(chain_prop)[i]=paste("Test_",colnames(chain_prop)[i],sep="")}
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    #Output chain to file every 10 iterations; start new file every 10,000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 0){fileIndex  = iter/10000}

      filename=paste(filename_prefix,fileIndex,".csv",sep="")
      lines=min((fileIndex * 10000+1),iter):iter
      cat("\nIteration ",iter,sep="")
      data_out<-cbind(posterior_current,posterior_prop,exp(chain),flag_accept,exp(chain_prop),chain_cov_all)[lines,]
      write.csv(data_out,filename,row.names=FALSE)
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }

    #Next iteration in chain
    out = MCMC_step(log_params,input_data,obs_sero_data,obs_case_data,chain_cov,adapt,like_current,
                     const_list)
  }

  #Get final parameter values
  param_out=exp(out$log_params)
  names(param_out)=names(log_params_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title MCMC_step
#'
#' @description Single MCMC step - one or more regions
#'
#' @details This function runs a single step in a Markov chain set up using the function mcmc(). It proposes a
#' set of parameters using the param_prop_setup() function, calculates the likelihood of observing the observed data
#' based on that proposed parameter set, accepts or rejects the proposed parameter set based on the calculated
#' likelihood and existing chain information, then returns the next line of information for the chain to mcmc().
#'
#' @param log_params Varied parameters (log values of actual parameters)
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param chain_cov = Chain covariance
#' @param adapt = 0/1 flag indicating which type of calculation to use for proposition value
#' @param like_current = Current accepted likelihood value
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc() (type,log_params_min,log_params_max,
#'   mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death,
#'   m_FOI_Brazil)
#'
#' @export
#'
MCMC_step <- function(log_params=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                       chain_cov=1,adapt=0,like_current=-Inf,const_list=list()) {

  #Propose new parameter values
  log_params_prop=param_prop_setup(log_params,chain_cov,adapt)

  #Calculate likelihood using single_like_calc function
  like_prop=single_like_calc(log_params_prop,input_data,obs_sero_data,obs_case_data,const_list)

  if(is.finite(like_prop)==FALSE) {
    p_accept = -Inf
  } else {
    p_accept= like_prop - like_current
    if(is.na(p_accept) ){ p_accept = -Inf}
  }

  ## accept/reject step:
  tmp = runif(1)
  if(tmp<min(exp(p_accept),1)) { # accept:
    log_params = log_params_prop
    like_current = like_prop
    accept = 1
  } else { # reject:
    accept = 0
  }

  return(list(log_params=log_params,log_params_prop=log_params_prop,like_current=like_current,like_prop=like_prop,
              accept=accept))
}
#-------------------------------------------------------------------------------
#' @title single_like_calc
#'
#' @description Function which calculates and outputs likelihood of observing simulated data
#'
#' @details This function calculates the total likelihood of observing a set of observations (across multiple
#' regions and data types) for a given proposed parameter set.
#'
#' @param log_params_prop Proposed values of varied parameters (natural logarithm of actual parameters)
#' @param input_data List of population and vaccination data for multiple regions (created using data input
#'   creation code and usually loaded from RDS file), with cross-reference tables added using input_data_process2
#'   in MCMC2
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc() (type,log_params_min,log_params_max,
#'   mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death,
#'   m_FOI_Brazil)
#'
#' @export
#'
single_like_calc <- function(log_params_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                              const_list=list()) {

  regions=input_data$region_labels
  n_regions=length(regions)

  #Get vaccine efficacy and calculate associated prior
  if(is.numeric(const_list$vaccine_efficacy)==FALSE){
    vaccine_efficacy=exp(log_params_prop[names(log_params_prop)=="vaccine_efficacy"])
    prior_vacc=log(dtrunc(vaccine_efficacy,"norm",a=0,b=1,mean=0.975,sd=0.05))
  } else {
    vaccine_efficacy=const_list$vaccine_efficacy
    prior_vacc=0
  }

  #Get reporting probabilities and check they are within specified bounds
  prior_report=0
  if(is.numeric(const_list$p_rep_severe)==FALSE){
    p_rep_severe=as.numeric(exp(log_params_prop[names(log_params_prop)=="p_rep_severe"]))
    if(p_rep_severe<exp(const_list$log_params_min[names(const_list$log_params_min)=="p_rep_severe"])){
      prior_report=-Inf}
    if(p_rep_severe>exp(const_list$log_params_max[names(const_list$log_params_max)=="p_rep_severe"])){
      prior_report=-Inf}
  } else {
    p_rep_severe=const_list$p_rep_severe
  }
  if(is.numeric(const_list$p_rep_death)==FALSE){
    p_rep_death=as.numeric(exp(log_params_prop[names(log_params_prop)=="p_rep_death"]))
    if(p_rep_death<exp(const_list$log_params_min[names(const_list$log_params_min)=="p_rep_death"])){
      prior_report=-Inf}
    if(p_rep_death>exp(const_list$log_params_max[names(const_list$log_params_max)=="p_rep_death"])){
      prior_report=-Inf}
  } else {
    p_rep_death=const_list$p_rep_death
  }

  #Get value of Brazil multiplier and check bounds
  if(is.numeric(const_list$m_FOI_Brazil)==FALSE){
    m_FOI_Brazil=as.numeric(exp(log_params_prop[names(log_params_prop)=="m_FOI_Brazil"]))
    if(m_FOI_Brazil<exp(const_list$log_params_min[names(const_list$log_params_min)=="m_FOI_Brazil"])){
      prior_report=-Inf}
    if(m_FOI_Brazil>exp(const_list$log_params_max[names(const_list$log_params_max)=="m_FOI_Brazil"])){
      prior_report=-Inf}
  } else {
    m_FOI_Brazil=const_list$m_FOI_Brazil
  }

  #Get FOI and R0 values and calculate associated prior
  FOI_R0_data=mcmc_FOI_R0_setup(const_list$type,const_list$prior_type,regions,log_params_prop,
                                const_list$enviro_data,const_list$R0_fixed_values,const_list$log_params_min,
                                const_list$log_params_max)
  FOI_values=FOI_R0_data$FOI_values
  for(n_region in 1:n_regions){
    if(substr(regions[n_region],1,3)=="BRA"){FOI_values[n_region]=FOI_values[n_region]*m_FOI_Brazil}
  }
  R0_values=FOI_R0_data$R0_values
  prior_prop=prior_vacc+prior_report+FOI_R0_data$prior+sum(dnorm(log(c(p_rep_severe,p_rep_death)),
                                                                 mean = 0,sd = 30,log = TRUE))

  if(is.null(obs_sero_data)){sero_like_values=NA}
  if(is.null(obs_case_data)){cases_like_values=deaths_like_values=NA}

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset(input_data,FOI_values,R0_values,obs_sero_data,obs_case_data,
                                vaccine_efficacy,p_rep_severe,p_rep_death,const_list$mode_start,const_list$dt)

    #Likelihood of observing serological data
    if(is.null(obs_sero_data)==FALSE){
      sero_like_values=lgamma(obs_sero_data$samples+1)-lgamma(obs_sero_data$positives+1)-
        lgamma(obs_sero_data$samples-obs_sero_data$positives+1)+
        obs_sero_data$positives*log(dataset$model_sero_values)+
        (obs_sero_data$samples-obs_sero_data$positives)*log(1.0-dataset$model_sero_values)
    }
    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data)==FALSE){
      model_case_values=dataset$model_case_values
      model_death_values=dataset$model_death_values
      for(i in 1:length(model_case_values)){
        model_case_values[i]=max(model_case_values[i],0.1)
        model_death_values[i]=max(model_death_values[i],0.1)
      }
      cases_like_values=dnbinom(x=obs_case_data$cases,mu=model_case_values,
                                size=rep(1,length(obs_case_data$cases)),log=TRUE)
      deaths_like_values=dnbinom(x=obs_case_data$deaths,mu=model_death_values,
                                 size=rep(1,length(obs_case_data$deaths)),log=TRUE)
    }

    likelihood=prior_prop+mean(c(sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
                     sum(deaths_like_values,na.rm=TRUE)),na.rm=TRUE)


  } else {likelihood=-Inf}

  return(likelihood)
}
