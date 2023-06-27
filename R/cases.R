# R file for functions relating to annual case and death data in YEP package
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
