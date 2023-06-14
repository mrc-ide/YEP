# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YEP package
#-------------------------------------------------------------------------------
#' @title MCMC2
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param log_params_ini Initial values of parameters to be varied. These should always be the log() values of the
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
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
#'   enviro"); set to NULL if not in use
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#' @param m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' '
#' @export
#'
MCMC2 <- function(log_params_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,filename_prefix="Chain",
                 Niter=1,type=NULL,log_params_min=c(),log_params_max=c(),mode_start=0,prior_type="zero",dt=1.0,
                 n_reps=1,enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,p_severe_inf=0.12,p_death_severe_inf=0.39,
                 p_rep_severe=NULL,p_rep_death=NULL,m_FOI_Brazil=1.0,deterministic=FALSE,mode_parallel="none",cluster=NULL){

  assert_that(is.logical(deterministic))

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
  consts=list(type=type,log_params_min=log_params_min,log_params_max=log_params_max,mode_start=mode_start,
              prior_type=prior_type,dt=dt,n_reps=n_reps,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
              vaccine_efficacy=vaccine_efficacy,p_severe_inf=p_severe_inf,p_death_severe_inf=p_death_severe_inf,
              p_rep_severe=p_rep_severe,p_rep_death=p_rep_death,m_FOI_Brazil=m_FOI_Brazil,deterministic=deterministic,
              mode_parallel=mode_parallel,cluster=cluster)


  #MCMC setup
  chain=chain_prop=posterior_current=posterior_prop=flag_accept=chain_cov_all=NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params=log_params_ini
  chain_cov=1
  adapt=0
  like_current=-Inf

  #Iterative fitting
  for (iter in 1:Niter){

    #cat("\n\tGenerating new parameter values")
    #Propose new parameter values
    log_params_prop=param_prop_setup(log_params,chain_cov,adapt)

    #cat("\n\tCalculating likelihood")
    #Calculate likelihood using single_like_calc function
    like_prop=single_like_calc(log_params_prop,input_data,obs_sero_data,obs_case_data,consts)
    #cat("\n\tLikelihood calculated")

    if(is.finite(like_prop)==FALSE) {
      p_accept = -Inf
    } else {
      p_accept = like_prop - like_current
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

    #save current step
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
      if(file.exists(filename)==FALSE){file.create(filename)}
      lines=min((fileIndex * 10000+1),iter):iter
      #cat("\nIteration ",iter,sep="")
      data_out<-cbind(posterior_current,posterior_prop,exp(chain),flag_accept,exp(chain_prop),chain_cov_all)[lines,]
      if(fileAccess(filename,2)==0){write.csv(data_out,filename,row.names=FALSE)}
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
  }

  #Get final parameter values
  param_out=exp(log_params)
  names(param_out)=names(log_params_ini)

  return(param_out)
}
