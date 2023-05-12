# R file for general functions in YEP package
#------------------------------------------------
#Global variables
p_severe_inf=0.12 #Probability that an infection is severe
p_death_severe_inf=0.39 #Probability that a severe infection becomes fatal
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YEP, .registration = TRUE
#' @importFrom assertthat assert_that
#' @import dde
#' @importFrom graphics axis matplot par
#' @importFrom mvtnorm rmvnorm
#' @import odin.dust
#' @import parallel
#' @importFrom R.utils fileAccess
#' @importFrom stats cov dexp dnbinom dnorm rbinom runif
#' @importFrom tgp lhs
#' @importFrom truncdist dtrunc
#' @importFrom utils write.csv
#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YEP", libpath)
}
#-------------------------------------------------------------------------------
#' @title Model_Run
#'
#' @description Run SEIRV model for single region (Model_Run_Multi_Region can be used to run multiple regions in parallel)
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output: "full" = SEIRVC + FOI for all steps, "case" = annual new infections (C),
#'   "sero" = annual SEIRV, "case+sero" = annual SEIRVC
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                      start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0,
                      vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  x <- SEIRV_Model$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                            vaccine_efficacy,start_SEIRV,dt),
                       time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  if(output_type=="full"){
    dimensions=c(N_age,n_particles,t_pts_out)
    output_data=list(day=x_res[1,1,],year=x_res[2,1,])
    output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
    output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
    output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
    output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
    output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
    output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
    output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)
  } else {
    n_years=length(years_data)
    dimensions=c(N_age,n_particles,n_years)
    output_data=list(year=years_data)
    if(output_type=="case+sero" || output_type=="sero"){
      output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=dimensions)
      for(n_year in 1:n_years){
        pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
        for(n_p in 1:n_particles){
          output_data$S[,n_p,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p,pts])
          output_data$E[,n_p,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p,pts])
          output_data$I[,n_p,n_year]=rowMeans(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p,pts])
          output_data$R[,n_p,n_year]=rowMeans(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p,pts])
          output_data$V[,n_p,n_year]=rowMeans(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p,pts])
        }
      }
    }
    if(output_type=="case+sero" || output_type=="case"){
      output_data$C=array(0,dim=c(n_particles,n_years))
      for(n_year in 1:n_years){
        pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
        for(n_p in 1:n_particles){
          output_data$C[n_p,n_year]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pts])
        }
      }
    }
  }

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Parameter setup
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin SEIRV model.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
parameter_setup <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                            years_data=c(1941:1942),mode_start=0,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0){

  assert_that(length(pop_data[,1])>1) #TODO - msg
  assert_that(length(pop_data[1,])>1) #TODO - msg
  n_years=length(pop_data[,1])-1
  N_age=length(pop_data[1,])
  assert_that(length(vacc_data[,1])==n_years+1,msg="Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1,])==N_age,msg="No. age groups in population and vaccination data must match")
  assert_that(mode_start %in% c(0,1,2),msg="mode_start must have value 0, 1 or 2")
  assert_that(vaccine_efficacy<=1.0 && vaccine_efficacy>=0.0,msg="Vaccine efficacy must be between 0 and 1")
  if(mode_start==2){
    assert_that(is.null(start_SEIRV$S)==FALSE,msg="When mode_start=2, start_SEIRV data is required")
  }
  assert_that(years_data[1]>=year0,msg="First data year must be greater than or equal to year0")
  assert_that(max(years_data)+1-years_data[1]<=n_years,msg="Period of years_data must lie within population data")
  vacc_initial=vacc_data[1,]
  assert_that(dt %in% c(1,2.5,5),msg="dt must have value 1, 2.5 or 5 days (must have integer no. points/year)")
  inv_365=1.0/365.0

  P0=Cas0=Sus0=Exp0=Inf0=Rec0=Vac0=rep(0,N_age)
  dP1_all=dP2_all=vacc_rates=array(rep(0,N_age*n_years),dim=c(N_age,n_years))
  for(i in 1:N_age){
    P0[i]=max(1.0,pop_data[1,i]) #Set all population values to nonzero minimum to avoid NaN values
  }
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i,n_year]=max(1.0,pop_data[n_year+1,i])*inv_365
      dP2_all[i,n_year]=max(1.0,pop_data[n_year,i])*inv_365
      if(i==1){
        vacc_rates[i,n_year]=vacc_data[n_year+1,i]*inv_365
      } else {
        vacc_rates[i,n_year]=max(0.0,vacc_data[n_year+1,i]-vacc_data[n_year,i-1])*inv_365
      }
    }
  }

  if(mode_start==0){
    Sus0=P0*(1.0-vacc_initial)
  }
  if(mode_start==1)
  {
    if(R0>1.0){
      herd_immunity=1.0-(1.0/R0)
    } else {
      herd_immunity=0.0
    }
    for(i in 1:N_age){
      if(vacc_initial[i]<herd_immunity){
        Rec0[i]=P0[i]*(herd_immunity-vacc_initial[i])
        Sus0[i]=P0[i]*(1.0-herd_immunity)
      } else {
        Sus0[i]=P0[i]*(1.0-vacc_initial[i])
      }
    }
  }
  if(mode_start==2){
    Sus0=start_SEIRV$S
    Exp0=start_SEIRV$E
    Inf0=start_SEIRV$I
    Rec0=start_SEIRV$R
    Vac0=start_SEIRV$V
    Cas0=rep(0,N_age)
  } else {
    Vac0=P0*vacc_initial
  }

  return(list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_annual=vacc_rates,
              Cas0=Cas0,Exp0=Exp0,Inf0=Inf0,N_age=N_age,Rec0=Rec0,Sus0=Sus0,Vac0=Vac0,
              dP1_all=dP1_all,dP2_all=dP2_all,n_years=n_years,year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt,
              t_incubation=t_incubation,t_latent=t_latent,t_infectious=t_infectious))
}
#-------------------------------------------------------------------------------
#' @title Generate_Dataset
#'
#' @description Generate serological and/or annual case/death data - TODO: add option(s) for parallelism
#'
#' @details This function is used to generate serological, annual case/death and/or annual outbreak risk data based
#'   on observed or dummy data sets; it is normally used by single_like_calc() and data_match_single() functions
#'
#' @param input_data List of population and vaccination data for multiple regions
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' '
#' @export
#'
Generate_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),obs_sero_data = NULL,obs_case_data = NULL,
                             vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,mode_start = 1,
                             start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE, mode_parallel = "none",
                             cluster = NULL){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                        " (see [TBA] )"))
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE),
              msg="Need at least one of obs_sero_data or obs_case_data")
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  if(is.null(obs_case_data)==FALSE){
    assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,
                msg="Severe infection reporting probability must be between 0 and 1")
    assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,
                msg="Fatal infection reporting probability must be between 0 and 1")
  }
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  if(is.null(input_data$flag_sero)){input_data=input_data_process(input_data,obs_sero_data,obs_case_data)}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)){model_sero_data=NULL} else {
    blank=rep(0,nrow(obs_sero_data))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  }
  if(is.null(obs_case_data)){model_case_values=model_death_values=NA} else {
    model_case_values=model_death_values=rep(0,nrow(obs_case_data))
  }

  #Set up vector of output types to get from model if needed
  if(mode_parallel %in% c("none","clusterMap")){
    output_types=rep(NA,n_regions)
    for(n_region in 1:n_regions){
      if(input_data$flag_case[n_region]==1){
        if(input_data$flag_sero[n_region]==1){output_types[n_region]="full"} else{output_types[n_region]="case"}
      } else {output_types[n_region]="sero"}
    }
  }

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(input_data$year_data_begin):max(input_data$year_end))
    if(any(input_data$flag_sero==1)){if(any(input_data$flag_case==1)){output_type="case+sero"} else{
      output_type="sero"}} else {output_type="case"}
    model_output_all=Model_Run_Multi_Region(FOI_spillover = FOI_values,R0 = R0_values,
                                            vacc_data = input_data$vacc_data,pop_data = input_data$pop_data,
                                            years_data = years_data_all,start_SEIRV=start_SEIRV,output_type = output_type,
                                            year0 = input_data$years_labels[1],mode_start = mode_start,
                                            vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                            n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV, output_type = output_types,
                                MoreArgs=list(year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(input_data$year_data_begin[n_region]:input_data$year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = output_types[n_region],
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
      #cat("\n\t\tFinished modelling region ",n_region)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile case data if needed
    if(input_data$flag_case[n_region]==1){
      case_line_list=input_data$case_line_list[[n_region]]
      years_outbreak=obs_case_data$year[case_line_list]
      n_years_outbreak=length(case_line_list)

      rep_cases=rep_deaths=rep(0,n_years_outbreak)
      for(n_year in 1:n_years_outbreak){ #TODO - Add functionality for n_reps>1
        pts=c(1:t_pts)[model_output$year==years_outbreak[n_year]]
        infs=model_output$C[pts]
        if(deterministic){
          severe_infs=floor(infs)*p_severe_inf
          deaths=severe_infs*p_death_severe_inf
          rep_deaths[n_year]=round(deaths*p_rep_death)
          rep_cases[n_year]=rep_deaths[n_year]+round((severe_infs-deaths)*p_rep_severe)

        } else {
          severe_infs=rbinom(1,floor(infs),p_severe_inf)
          deaths=rbinom(1,severe_infs,p_death_severe_inf)
          rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
          rep_cases[n_year]=rep_deaths[n_year]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
        }
      }

      model_case_values[case_line_list]=model_case_values[case_line_list]+rep_cases
      model_death_values[case_line_list]=model_death_values[case_line_list]+rep_deaths
    }

    #Compile seroprevalence data if necessary
    if(input_data$flag_sero[n_region]==1){
      sero_line_list=input_data$sero_line_list[[n_region]]
      sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_output)
      model_sero_data$samples[sero_line_list]=model_sero_data$samples[sero_line_list]+sero_results$samples
      model_sero_data$positives[sero_line_list]=model_sero_data$positives[sero_line_list]+sero_results$positives
    }
    model_output<-NULL
  }
  model_output_all<-NULL

  if(any(input_data$flag_sero>0)){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(any(input_data$flag_case>0)){
    model_case_values=model_case_values/n_reps
    model_death_values=model_death_values/n_reps
  }

  return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
              model_death_values=model_death_values))
}
#-------------------------------------------------------------------------------
#' @title create_param_labels
#'
#' @description Apply names to the parameters in a set used for data matching and parameter fitting
#'
#' @details Takes in input list and environmental data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#' code and usually loaded from RDS file)
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param extra_params Vector of strings listing parameters besides ones determining FOI/R0 (may include vaccine
#'   efficacy and/or infection/death reporting probabilities)
#'
#' @export
#'
create_param_labels <- function(type="FOI",input_data=list(),enviro_data=NULL,extra_params=c("vacc_eff")){

  assert_that(type %in% c("FOI","FOI+R0","FOI enviro","FOI+R0 enviro"))
  assert_that(input_data_check(input_data),msg="Input data must be in standard format")

  n_extra=length(extra_params)

  if(type %in% c("FOI","FOI+R0")){
    regions=input_data$region_labels
    n_regions=length(regions)
    if(type=="FOI"){n_params=n_regions+n_extra} else {n_params=(2*n_regions)+n_extra}
    param_names=rep("",n_params)
    for(i in 1:n_regions){
      param_names[i]=paste("FOI_",regions[i],sep="")
      if(type=="FOI+R0"){param_names[i+n_regions]=paste("R0_",regions[i],sep="")}
    }
  } else {
    assert_that(is.data.frame(enviro_data)) #TODO - msg
    env_vars=colnames(enviro_data)[c(2:ncol(enviro_data))]
    n_env_vars=length(env_vars)
    if(type=="FOI enviro"){n_params=n_env_vars+n_extra} else {n_params=(2*n_env_vars)+n_extra}
    param_names=rep("",n_params)
    for(i in 1:n_env_vars){
      param_names[i]=paste("FOI_",env_vars[i],sep="")
      if(type=="FOI+R0 enviro"){param_names[i+n_env_vars]=paste("R0_",env_vars[i],sep="")}
    }
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params]=extra_params}

  return(param_names)
}
#-------------------------------------------------------------------------------
#' @title param_calc_enviro
#'
#' @description Parameter calculation from environmental covariates
#'
#' @details Takes in set of coefficients of environmental covariates and covariate values and calculates values of
#'   spillover force of infection and reproduction number.
#'
#' @param enviro_coeffs Values of environmental coefficients
#' @param enviro_covar_values Values of environmental covariates
#' '
#' @export
#'
param_calc_enviro <- function(enviro_coeffs=c(),enviro_covar_values=c()){

  assert_that(all(enviro_coeffs>=0),msg="All environmental coefficients must have positive values")
  n_env_vars=length(enviro_covar_values)
  assert_that(length(enviro_coeffs) %in% c(n_env_vars,2*n_env_vars),msg="Wrong number of environmental coefficients")

  output=list(FOI=NA,R0=NA)
  output$FOI=sum(enviro_coeffs[c(1:n_env_vars)]*enviro_covar_values)
  if(length(enviro_coeffs)==2*n_env_vars){
    output$R0=sum(enviro_coeffs[c(1:n_env_vars)+n_env_vars]*enviro_covar_values)
    }

  return(output)
}
