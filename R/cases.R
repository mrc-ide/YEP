# R file for functions relating to annual case and death data in YEP package
#-------------------------------------------------------------------------------
#' @title case_data_calculate
#'
#' @description Calculate reported case data from SEIRV model output
#'
#' @details [TBA]
#'
#' @param model_data
#' @param n_p Particle to select from model_data
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' '
#' @export
#'
case_data_calculate <- function(model_data = list(), n_p = 1, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                p_rep_severe = 1.0, p_rep_death = 1.0){
  #TODO - Add assert_that functions
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
