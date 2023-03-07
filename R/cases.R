# R file for functions relating to annual case and death data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title case_data_generate
#'
#' @description Take in single set of population data and model parameters, output infection data only
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period and outputs infection numbers at each
#' time increment only; optimized for running a large number of repetitions
#'
#' @param FOI_spillover = Force of infection due to spillover from sylvatic reservoir
#' @param R0 = Basic reproduction number for urban spread of infection
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
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' '
#' @export
#'
case_data_generate <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,years_data=c(),
                               mode_start=0,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0) {

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,start_SEIRV,dt)

  x <- SEIRV_Model$new(FOI_spillover=pars$FOI_spillover,R0=pars$R0,vacc_rate_annual=pars$vacc_rate_annual,
                       Cas0=pars$Cas0,Exp0=pars$Exp0,Inf0=pars$Inf0,N_age=pars$N_age,Rec0=pars$Rec0,Sus0=pars$Sus0,
                       Vac0=pars$Vac0,dP1_all=pars$dP1_all,dP2_all=pars$dP2_all,n_years=pars$n_years,year0=pars$year0,
                       vaccine_efficacy=pars$vaccine_efficacy,dt=pars$dt)

  n_nv=5 #Number of non-vector outputs
  N_age=length(pop_data[1,])
  t_pts=c(1:(((max(years_data)+1)-year0)*(365/dt)))
  n_data_pts=(6*N_age)+n_nv
  n_steps=length(t_pts)
  step0=(years_data[1]-year0)*(365/dt)
  steps=n_steps-step0
  results_data=list(year=sort(rep(c(years_data[1]:max(years_data)),(365/dt))),C=rep(0,steps))
  pts_select=c(((5*N_age)+1+n_nv):((6*N_age)+n_nv))

  x_res <- x$run(n_steps)
  for(t in c((step0+1):n_steps)){
    results_data$C[t-step0]=sum(x_res[t,pts_select])
  }

  return(results_data)
}
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
#' @title deaths_compare
#'
#' @description Compare modelled and observed deaths using negative binomial
#'
#' @details Compares modelled data on fatal cases per year and compared with observed data, calculating logarithmic
#' likelihood of observing the latter given the former, using a negative binomial formula.
#'
#' @param model_data Modelled data in data frame format (list of years and number of reported deaths in each)
#' @param obs_data Data frame containing observed death data
#' '
#' @export
#'
deaths_compare <- function(model_data=list(),obs_data=list()){

  assert_that(is.null(model_data$rep_deaths)==FALSE,msg="Modelled data must include reported deaths")
  assert_that(is.null(obs_data$deaths)==FALSE,msg="Observed data must include numbers of deaths")
  assert_that(length(model_data$rep_deaths)==length(obs_data$deaths),
              msg="Numbers of entries in observed and modelled data must match")
  model_data$rep_deaths[model_data$rep_deaths==0]=0.1

  like_values=dnbinom(x=obs_data$deaths,mu=model_data$rep_deaths,size=rep(1,length(obs_data$deaths)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
#' @title cases_compare
#'
#' @description Compare modelled and observed severe cases using negative binomial
#'
#' @details Compares modelled data on severe cases per year and compared with observed data, calculating logarithmic
#' likelihood of observing the latter given the former, using a negative binomial formula.
#'
#' @param model_data Modelled data in data frame format (list of years and number of reported severe cases in each)
#' @param obs_data Data frame containing annual observed case data
#' '
#' @export
#'
cases_compare <- function(model_data=list(),obs_data=list()){

  assert_that(is.null(model_data$rep_cases)==FALSE,msg="Modelled data must include reported cases")
  assert_that(is.null(obs_data$cases)==FALSE,msg="Observed data must include reported cases")
  assert_that(length(model_data$rep_cases)==length(obs_data$cases),
              msg="Numbers of entries in modelled and observed data must match")
  model_data$rep_cases[model_data$rep_cases==0]=0.1

  like_values=dnbinom(x=obs_data$cases,mu=model_data$rep_cases,size=rep(1,length(obs_data$cases)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
