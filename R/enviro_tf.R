#-------------------------------------------------------------------------------
#' @title briere
#'
#' @description Briere function for temperature transformation
#'
#' @details TBA
#'
#' @param temp_values TBA
#' @param T0 TBA
#' @param Tm TBA
#' @param c TBA
#' '
#' @export
#'
briere = function(temp_values = c(25.0), T0 = 2.248952, Tm = 40.13383, c = 0.000271964){
  out = c*temp_values*(temp_values-T0)*(Tm-temp_values)^0.5
  return(out)
}
#-------------------------------------------------------------------------------
#' @title quad
#'
#' @description Quadratic function for temperature transformation
#'
#' @details TBA
#'
#' @param temp_values TBA
#' @param T0 TBA
#' @param Tm TBA
#' @param c TBA
#' '
#' @export
#'
quad = function(temp_values = c(25.0), T0 = 12.71508, Tm = 38.04809469, c = -0.757869){
  out = -c*(temp_values-T0)*(Tm -temp_values)
  return(out)
}
#-------------------------------------------------------------------------------
#' @title temp_tf
#'
#' @description  Function for temperature transformation
#'
#' @details TBA
#'
#' @param temp_values TBA
#' @param a_T0 TBA
#' @param a_Tm TBA
#' @param a_c TBA
#' @param mu_T0 TBA
#' @param mu_Tm TBA
#' @param mu_c TBA
#' @param PDR_T0 TBA
#' @param PDR_Tm TBA
#' @param PDR_c TBA
#' '
#' @export
#'
temp_tf <- function(temp_values = c(25.0),
                    a_T0 = 2.248952, a_Tm = 40.13383, a_c = 0.000271964,
                    mu_T0 = 12.71508, mu_Tm = 38.04809469, mu_c = -0.757869,
                    PDR_T0 = 17.33263, PDR_Tm = 42.19592, PDR_c = 0.000135891){


  #bite rate
  a = briere(temp_values, a_T0, a_Tm, a_c) # new estimate medians

  # mu = m mortality
  lf = quad(temp_values, mu_T0, mu_Tm, mu_c)
  mu = 1/ ifelse( lf<=0, 1, lf ) #guard against negatives and zeros

  # PDR = parasite development rate = 1/EIP
  PDR = briere(temp_values, PDR_T0, PDR_Tm, PDR_c) #
  PDR[PDR<0]=0

  a[is.na(a)] = PDR[is.na(PDR)] = 0

  # suitability
  temp_tf_values = (a^2 * exp(-mu / PDR) ) / mu

  return(temp_tf_values)
}
#-------------------------------------------------------------------------------
#' @title precip_tf
#'
#' @description Function for precipitation transformation
#'
#' @details TBA
#'
#' @param precip_values TBA
#' @param a_ptf TBA
#' '
#' @export
#'
precip_tf <- function(precip_values = c(0), a_ptf = 0){
  p_lim=0.5/a_ptf
  precip_values[precip_values>p_lim]=p_lim
  precip_tf_values = precip_values*(1.0 - (a_ptf*precip_values))

  return(precip_tf_values)
}
#-------------------------------------------------------------------------------
#' @title epi_param_calc2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param pars_fixed TBA
#' @param env_covar_values TBA
#' @param log_FOI_coeffs TBA
#' @param log_R0_coeffs TBA
#' @param vars_extra TBA
#'
#' @export
#'
epi_param_calc2 <- function(pars_fixed = list(), env_covar_values = list(), log_FOI_coeffs = c(),
                            log_R0_coeffs = c(), vars_extra = list()){

  #TODO - adapt for multiple parameter sets?

  if("m_FOI_BRA" %in% names(vars_extra)){flag_BRA = 2} else {
    if(is.null(pars_fixed$m_FOI_BRA) == FALSE){ flag_BRA = 1 } else {flag_BRA = 0}
  }
  n_regions = pars_fixed$n_r
  assert_that(n_regions==dim(env_covar_values)[2])
  dim_t = dim(env_covar_values)[3]
  time_inc = pars_fixed$time_inc
  pts_year = 365.0/time_inc
  n_years = pars_fixed$n_years
  n_t_pts = n_years*pts_year
  inv_365 = 1.0/365.0

  a_T0 = a_Tm = a_c = mu_T0 = mu_Tm = mu_c = PDR_T0 = PDR_Tm = PDR_c = log_a_ptf = 0
  #Temperature; TODO - add functionality for vector i_ttf (multiple temperature covariates)
  if(is.null(pars_fixed$i_ttf)==FALSE){
    for(name in extra_param_names_ttf){
      if(is.null(pars_fixed[[name]])){assign(name,vars_extra[[name]])}else{assign(name,pars_fixed[[name]])}
    }
    env_covar_values[pars_fixed$i_ttf,,]=temp_tf(temp_values=array(env_covar_values[pars_fixed$i_ttf,,],
                                                                   dim=c(n_regions,dim_t)),
                                                 a_T0, a_Tm, a_c,mu_T0, mu_Tm, mu_c, PDR_T0, PDR_Tm,PDR_c)
  }
  #Precipitation; TODO - ditto
  if(is.null(pars_fixed$i_ptf)==FALSE){
    if(is.null(pars_fixed$log_a_ptf)){log_a_ptf=vars_extra[["log_a_ptf"]]}else{log_a_ptf=pars_fixed$log_a_ptf}
    env_covar_values[pars_fixed$i_ptf,,] = precip_tf(precip_values=array(env_covar_values[pars_fixed$i_ptf,,],
                                                                         dim=c(n_regions,dim_t)),
                                                     a_ptf = exp(log_a_ptf))
  }

  FOI_spillover = array(colSums(as.numeric(exp(log_FOI_coeffs))*env_covar_values),dim=c(n_regions,dim_t))
  R0 = array(colSums(as.numeric(exp(log_R0_coeffs))*env_covar_values),dim=c(n_regions,dim_t))
  if(flag_BRA>0){
    if(flag_BRA==1){m=pars_fixed$m_FOI_BRA}else{m=vars_extra$m_FOI_BRA}
    FOI_spillover[pars_fixed$ref_BRA,]=FOI_spillover[pars_fixed$ref_BRA,]*m
  }

  return(list(FOI_spillover = FOI_spillover,R0 = R0))
}
