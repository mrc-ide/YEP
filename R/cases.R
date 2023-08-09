# R file for functions relating to annual case and death data in YEP package
#-------------------------------------------------------------------------------
#' @title case_data_calculate
#'
#' @description Calculate reported case data from SEIRV model output
#'
#' @details [TBA]
#'
#' @param model_data SEIRV output of Model_Run and similar functions
#' @param n_p Particle to select from model_data
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param output_type Type of output to produce:
#'   "annual" - Total reported cases and reported deaths by year
#'   "pts" - Total reported cases and reported deaths for every time point, summed over age groups
#'   "full" - Reported cases and reported deaths for every time point and age group
#' @param deterministic Indicates whether to calculate results deterministically (TRUE) or stochastically (FALSE)
#' '
#' @export
#'
case_data_calculate <- function(model_data = list(), n_p = 1, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                p_rep_severe = 1.0, p_rep_death = 1.0, output_type = "annual", deterministic = FALSE){
  #TODO - Add assert_that functions
  assert_that(length(dim(model_data$C)) %in% c(2,3))

  if(length(dim(model_data$C))==2){
    input_type="combined"
    assert_that(output_type %in% c("annual","pts"))
  } else {
    input_type="age_split"
    assert_that(output_type %in% c("annual","pts","full"))
  }

  if(output_type=="annual"){
    years=unique(model_data$year)
    n_years=length(years)
    case_data=data.frame(rep_cases=rep(NA,n_years),rep_deaths=rep(NA,n_years))
    for(n_year in 1:n_years){
      if(input_type=="age_split"){
        infs=sum(model_data$C[,n_p,model_data$year==years[n_year]])
      }else{
        infs=sum(model_data$C[n_p,model_data$year==years[n_year]])
      }
      severe_infs=rbinom(1,floor(infs),p_severe_inf)
      deaths=rbinom(1,severe_infs,p_death_severe_inf)
      case_data$rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
      case_data$rep_cases[n_year]=case_data$rep_deaths[n_year] + rbinom(1,severe_infs-deaths,p_rep_severe)
    }
  } else {
    if(output_type=="pts"){
      if(input_type=="age_split"){n_pts=dim(model_data$C)[3]}else{n_pts=dim(model_data$C)[2]}
      case_data=data.frame(rep_cases=rep(NA,n_pts),rep_deaths=rep(NA,n_pts))
      for(n_pt in 1:n_pts){
        if(input_type=="age_split"){infs=sum(model_data$C[,n_p,n_pt])}else{infs=sum(model_data$C[n_p,n_pt])}
        severe_infs=rbinom(1,floor(infs),p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        case_data$rep_deaths[n_pt]=rbinom(1,deaths,p_rep_death)
        case_data$rep_cases[n_pt]=case_data$rep_deaths[n_pt] + rbinom(1,severe_infs-deaths,p_rep_severe)
      }
    } else {
      N_age=dim(model_data$C)[1]
      case_data=list(rep_cases=array(NA,dim=c(N_age,n_pts)),rep_deaths=array(NA,dim=c(N_age,n_pts)))
      for(n_pt in 1:n_pts){
        if(input_type=="age_split"){infs=sum(model_data$C[,n_p,n_pt])}else{infs=sum(model_data$C[n_p,n_pt])}
        severe_infs=rbinom(1,floor(infs),p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        case_data$rep_deaths[,n_pt]=rbinom(1,deaths,p_rep_death)
        case_data$rep_cases[,n_pt]=case_data$rep_deaths[n_pt] + rbinom(1,severe_infs-deaths,p_rep_severe)
      }
    }
  }

  return(case_data)
}
#-------------------------------------------------------------------------------
#' @title case_data_calculate_multi
#'
#' @description Calculate reported case data from SEIRV model output across multiple repetitions
#'
#' @details [TBA]
#'
#' @param model_data SEIRV output of Model_Run and similar functions
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param output_type Type of output to produce:
#'   "annual" - Total reported cases and reported deaths by year
#'   "pts" - Total reported cases and reported deaths for every time point, summed over age groups
#'   "full" - Reported cases and reported deaths for every time point and age group
#' @param deterministic Indicates whether to calculate results deterministically (TRUE) or stochastically (FALSE)
#' '
#' @export
#'
case_data_calculate_multi <- function(model_data = list(), p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                      p_rep_severe = 1.0, p_rep_death = 1.0, output_type = "annual", deterministic = FALSE){
  #TODO - Add assert_that functions
  assert_that(length(dim(model_data$C)) %in% c(2,3))
  if(length(dim(model_data$C))==2){
    input_type="combined"
    assert_that(output_type %in% c("annual","pts"))
    n_reps=dim(model_data$C)[1]
  } else {
    input_type="age_split"
    assert_that(output_type %in% c("annual","pts","full"))
    n_reps=dim(model_data$C)[2]
  }

  case_data_multi=list()
  for(n_rep in 1:n_reps){
    case_data_multi[[n_rep]] = case_data_calculate(model_data, n_p = n_rep, p_severe_inf, p_death_severe_inf,
                                           p_rep_severe, p_rep_death,output_type, deterministic)
  }

  return(case_data_multi)
}
#-------------------------------------------------------------------------------
#' @title case_data_compare
#'
#' @description Compare modelled and observed case or death data using negative binomial function
#'
#' @details Compares modelled data (from dataset generation functions) on reported cases or deaths
#' per year with observed data, calculating logarithmiclikelihood of observing the latter given the former,
#' using a negative binomial formula.
#'
#' @param model_values Modelled reported case or death values
#' @param obs_values Observed or template case or death values
#' '
#' @export
#'
case_data_compare <- function(model_values,obs_values){

  #TODO - Add assert_that functions

  model_values[model_values<0.1]=0.1 #Eliminate values of 0 to avoid any -Inf values in like_values

  like_values=dnbinom(x=obs_values,mu=model_values,size=rep(1,length(obs_values)),log=TRUE)

  return(like_values)
}
