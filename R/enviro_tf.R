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
