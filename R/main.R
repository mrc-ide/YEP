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
#' @import assertthat
#' @import dde
#' @import ggplot2
#' @import graphics
#' @import grDevices
#' @import maptree
#' @import sf
#' @import stats
#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YEP", libpath)
}
#-------------------------------------------------------------------------------
#' @title Model_Run
#'
#' @description Run SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param years_data Incremental vector of years for which to output SEIRV data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRVC input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
Model_Run <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,years_data=c(1941:1942),
                      mode_start=0,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0) {

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,
                       start_SEIRV,dt)

  x <- SEIRV_Model$new(FOI_spillover=pars$FOI_spillover,R0=pars$R0,vacc_rate_annual=pars$vacc_rate_annual,
                       Cas0=pars$Cas0,Exp0=pars$Exp0,Inf0=pars$Inf0,N_age=pars$N_age,Rec0=pars$Rec0,Sus0=pars$Sus0,
                       Vac0=pars$Vac0,dP1_all=pars$dP1_all,dP2_all=pars$dP2_all,n_years=pars$n_years,
                       year0=pars$year0,vaccine_efficacy=pars$vaccine_efficacy,dt=pars$dt)

  n_nv=5 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  t_pts_all=c(1:((max(years_data)+1-year0)*(365/dt))) #All output time points
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  n_steps=length(t_pts_all) #Total number of output time points
  step0=(years_data[1]-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data
  x_res <- x$run(n_steps)
  t_pts=c((step0+1):n_steps)

  return(list(day=x_res[t_pts,3],year=x_res[t_pts,4],FOI_total=x_res[t_pts,5],
              S=array(x_res[t_pts,c((1+n_nv):(N_age+n_nv))],dim=c(t_pts_out,N_age)),
              E=array(x_res[t_pts,c((N_age+1+n_nv):((2*N_age)+n_nv))],dim=c(t_pts_out,N_age)),
              I=array(x_res[t_pts,c(((2*N_age)+1+n_nv):((3*N_age)+n_nv))],dim=c(t_pts_out,N_age)),
              R=array(x_res[t_pts,c(((3*N_age)+1+n_nv):((4*N_age)+n_nv))],dim=c(t_pts_out,N_age)),
              V=array(x_res[t_pts,c(((4*N_age)+1+n_nv):((5*N_age)+n_nv))],dim=c(t_pts_out,N_age)),
              C=array(x_res[t_pts,c(((5*N_age)+1+n_nv):((6*N_age)+n_nv))],dim=c(t_pts_out,N_age))))
}
#-------------------------------------------------------------------------------
#' @title Parameter setup
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin.dust SEIRV model versions.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param years_data Incremental vector of years for which to output SEIRV data
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
  assert_that(is.integer(years_data) &&
                all(years_data[c(2:length(years_data))]-years_data[c(1:(length(years_data)-1))]),
              msg="years_data must be iterative vector of years")
  assert_that(max(years_data)+1-years_data[1]<=n_years,msg="Period of years_data must lie within population data")
  vacc_initial=vacc_data[1,]
  assert_that(dt %in% c(1,2.5,5),msg="dt must have value 1, 2.5 or 5 days (must have integer no. points/year")
  inv_365=1.0/365.0

  P0=Cas0=Sus0=Exp0=Inf0=Rec0=Vac0=rep(0,N_age)
  dP1_all=dP2_all=vacc_rates=array(rep(0,N_age*n_years),dim=c(N_age,n_years))
  for(i in 1:N_age){
    P0[i]=max(1.0,pop_data[1,i]) #Set all population values to a minimum of 1 to avoid NaN values appearing
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
              dP1_all=dP1_all,dP2_all=dP2_all,n_years=n_years,year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt))
}
#-------------------------------------------------------------------------------
#' @title Generate_Dataset
#'
#' @description Generate serological, annual case/death and/or annual outbreak risk data
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
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
Generate_Dataset <- function(input_data=list(),FOI_values=c(),R0_values=c(),
                             obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                             vaccine_efficacy=1.0,p_rep_severe=1.0,p_rep_death=1.0,mode_start=1,dt=1.0){

  assert_that(input_data_check(input_data),
              msg=paste("Input data must be in standard format",
                        " (see https://mrc-ide.github.io/YellowFeverDynamics/articles/CGuideAInputs.html )"))
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE,is.null(obs_outbreak_data)==FALSE),
              msg="Need at least one of obs_sero_data, obs_case_data or obs_outbreak_data")
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  if(is.null(obs_case_data)==FALSE){
    assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,
                msg="Severe infection reporting probability must be between 0 and 1")
    assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,
                msg="Fatal infection reporting probability must be between 0 and 1")
    }

  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")

  if(is.null(input_data$flag_sero)){
    input_data=input_data_process(input_data,obs_sero_data,obs_case_data)
  }

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)==FALSE){
    blank=rep(0,nrow(obs_sero_data))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  } else {model_sero_data=NULL}

  if(is.null(obs_case_data)==FALSE){
    model_case_values=model_death_values=rep(0,nrow(obs_case_data))
  } else {model_case_values=model_death_values=NA}

  #Model all regions and save relevant output data
  for(n_region in 1:n_regions){

    #Get information on which observed data types are available for considered region
    flag_sero=input_data$flag_sero[n_region]
    flag_case=input_data$flag_case[n_region]

    #Get input data on region
    vacc_data=input_data$vacc_data[n_region,,]
    pop_data=input_data$pop_data[n_region,,]
    year_end=input_data$year_end[n_region]
    year_data_begin=input_data$year_data_begin[n_region]

    #Run model
    model_output=Model_Run(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                           year0=min(input_data$years_labels),years_data=c(year_data_begin:year_end),
                           mode_start,vaccine_efficacy,dt=dt)

    #Compile case data if needed
    if(flag_case==1){
      case_line_list=input_data$case_line_list[[n_region]]
      years_outbreak=obs_case_data$year[case_line_list]
      n_years_outbreak=length(case_line_list)

      blank=rep(0,n_years_outbreak)
      annual_data=list(rep_cases=blank,rep_deaths=blank)
      for(n_year in 1:n_years_outbreak){
        year=years_outbreak[n_year]
        infs=sum(model_output$C[model_output$year==year,])
        severe_infs=rbinom(1,floor(infs),p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        annual_data$rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
        annual_data$rep_cases[n_year]=annual_data$rep_deaths[n_year]+rbinom(1,severe_infs-deaths,p_rep_severe)
      }

      model_case_values[case_line_list]=annual_data$rep_cases
      model_death_values[case_line_list]=annual_data$rep_deaths
    }

    #Compile seroprevalence data if necessary
    if(flag_sero==1){
      sero_line_list=input_data$sero_line_list[[n_region]]
      sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_output,i)
      model_sero_data$samples[sero_line_list]=sero_results$samples
      model_sero_data$positives[sero_line_list]=sero_results$positives
    }
    model_output<-NULL
  }

  if(flag_sero==1){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}

  return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
              model_death_values=model_death_values))
}
