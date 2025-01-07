# R file for general functions in YEP package
#-------------------------------------------------------------------------------
# Global variables
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#-------------------------------------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YEP, .registration = TRUE
#' @importFrom assertthat assert_that
#' @import dde
#' @import dust2
#' @importFrom graphics axis matplot par
#' @importFrom mvtnorm rmvnorm
#' @import odin2
#' @import parallel
#' @importFrom R.utils fileAccess
#' @importFrom stats cov dexp dnbinom dnorm nlm rbinom runif
#' @importFrom tgp lhs
#' @importFrom truncdist dtrunc
#' @importFrom utils write.csv
#-------------------------------------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YEP", libpath)
}
#-------------------------------------------------------------------------------
#' @title Model_Run
#'
#' @description Run SEIRV model for single region
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values. Run up to 20 particles (Model_Run_Many_Reps
#' can be used to run larger numbers of particles).
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param output_type Type of data to output: \cr
#'   "full" = SEIRVC + FOI for all steps and ages \cr
#'   "case" = annual total new infections (C) summed across all ages \cr
#'   "sero" = annual SEIRV \cr
#'   "case+sero" = annual SEIRVC, C summed across all ages \cr
#'   "case_alt" = annual total new infections not combined by age \cr
#'   "case_alt2" = total new infections combined by age for all steps
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, output_type = "full", mode_start = 0,
                      start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB  -  Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1] - year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  x <- dust_system_create(SEIRV_Model, pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                                                              vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time),
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = c(step_begin:step_end))

  if(output_type == "full"){
    dim = c(N_age, n_particles, t_pts_out)
    output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ], FOI_total = x_res[3, , ]/time_inc,
                       S = array(x_res[index$S, , ], dim), E = array(x_res[index$E, , ], dim), I = array(x_res[index$I, , ], dim),
                       R = array(x_res[index$R, , ], dim), V = array(x_res[index$V, , ], dim), C = array(x_res[index$C, , ], dim))
  } else {
    if(output_type == "case_alt2"){
      output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ])
      output_data$C = colSums(x_res[index$C, , ])
    }  else {
      n_years = length(years_data)
      output_data = list(year = years_data)
      if(output_type == "case+sero" || output_type == "sero"){
        output_data$V = output_data$R = output_data$I = output_data$E = output_data$S = array(0, dim = c(N_age, n_particles, n_years))
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$S[, n_p, n_year] = rowMeans(x_res[index$S, n_p, pts])
            output_data$E[, n_p, n_year] = rowMeans(x_res[index$E, n_p, pts])
            output_data$I[, n_p, n_year] = rowMeans(x_res[index$I, n_p, pts])
            output_data$R[, n_p, n_year] = rowMeans(x_res[index$R, n_p, pts])
            output_data$V[, n_p, n_year] = rowMeans(x_res[index$V, n_p, pts])
          }
        }
      }
      if(output_type == "case+sero" || output_type == "case"){
        output_data$C = array(0, dim = c(n_particles, n_years))
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){ output_data$C[n_p, n_year] = sum(x_res[index$C, n_p, pts]) }
        }
      }
      if(output_type == "case_alt"){
        output_data$C = array(0, dim = c(N_age, n_particles, n_years))
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){ output_data$C[, n_p, n_year] = rowSums(x_res[index$C, n_p, pts]) }
        }
      }
    }
  }
  x_res = NULL
  gc()

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Many_Reps
#'
#' @description Run SEIRV model for single region for large number of repetitions
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of repetitions and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values. Variation of Model_Run() used for
#' running a large number of repetitions (>20).
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param output_type Type of data to output: \cr
#'   "full" = SEIRVC + FOI for all steps and ages \cr
#'   "case" = annual total new infections (C) summed across all ages \cr
#'   "sero" = annual SEIRV \cr
#'   "case+sero" = annual SEIRVC, C summed across all ages \cr
#'   "case_alt" = annual total new infections not combined by age \cr
#'   "case_alt2" = total new infections combined by age for all steps
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param n_reps Number of repetitions (used to set number of particles and threads)
#' @param division Number of particles/threads to run in one go (up to 20)
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Many_Reps <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                                year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, output_type = "full", mode_start = 0,
                                start_SEIRV = list(), mode_time = 0, n_reps = 1, division = 10, deterministic = FALSE) {

  assert_that(division <= 20, msg = "Number of particles run at once must be 20 or less")
  n_particles0 = min(division, n_reps)
  n_threads = min(division, n_particles0)
  n_divs = ceiling(n_reps/division)
  if(n_divs == 1){
    n_particles_list = n_particles0
  } else {
    n_particles_list = c(rep(n_particles0, n_divs - 1), n_reps - (division*(n_divs - 1)))
  }

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1] - year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  if(output_type == "full"){
    dim = c(N_age, n_reps, t_pts_out)
    output_data = list(day = rep(NA, t_pts_out), year = rep(NA, t_pts_out), FOI_total = array(NA, c(n_reps, t_pts_out)),
                       S = array(NA, dim), E = array(NA, dim), I = array(NA, dim),
                       R = array(NA, dim), V = array(NA, dim), C = array(NA, dim))
  } else {
    if(output_type == "case_alt2"){
      output_data = list(day = rep(NA, t_pts_out), year = rep(NA, t_pts_out))
      output_data$C = array(0, dim = c(n_reps, t_pts_out))
    } else {
      n_years = length(years_data)
      output_data = list(year = years_data)
      if(output_type == "case+sero" || output_type == "sero"){
        output_data$V = output_data$R = output_data$I = output_data$E = output_data$S = array(0, dim = c(N_age, n_reps, n_years))
      }
      if(output_type == "case+sero" || output_type == "case"){ output_data$C = array(0, dim = c(n_reps, n_years)) }
      if(output_type == "case_alt"){ output_data$C = array(0, dim = c(N_age, n_reps, n_years)) }
    }
  }

  pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                         vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)
  for(div in 1:n_divs){
    n_particles = n_particles_list[div]

    x <- dust_system_create(SEIRV_Model, pars = pars, n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                            deterministic = FALSE, preserve_particle_dimension = TRUE)
    dust_system_set_state_initial(x)
    x_res <- dust_system_simulate(x, times = c(step_begin:step_end))

    if(div == 1){
      n_p0 = 0
      index = dust_unpack_index(x)
    } else{
      n_p0 = sum(n_particles_list[c(1:(div - 1))])
    }

    if(output_type == "full"){
      n_p_values = c(1:n_particles) + n_p0
      dim = c(N_age, n_particles, t_pts_out)
      output_data$day = x_res[1, 1, ]
      output_data$year = x_res[2, 1, ]
      output_data$FOI_total[n_p_values, ] = array(x_res[3, , ]/time_inc, dim = c(n_particles, t_pts_out))
      output_data$S[, n_p_values, ] = array(x_res[index$S, , ], dim)
      output_data$E[, n_p_values, ] = array(x_res[index$E, , ], dim)
      output_data$I[, n_p_values, ] = array(x_res[index$I, , ], dim)
      output_data$R[, n_p_values, ] = array(x_res[index$R, , ], dim)
      output_data$V[, n_p_values, ] = array(x_res[index$V, , ], dim)
      output_data$C[, n_p_values, ] = array(x_res[index$C, , ], dim)
    } else {
      if(output_type == "case+sero" || output_type == "sero"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p + n_p0
            output_data$S[, n_p2, n_year] = rowMeans(x_res[index$S, n_p, pts])
            output_data$E[, n_p2, n_year] = rowMeans(x_res[index$E, n_p, pts])
            output_data$I[, n_p2, n_year] = rowMeans(x_res[index$I, n_p, pts])
            output_data$R[, n_p2, n_year] = rowMeans(x_res[index$R, n_p, pts])
            output_data$V[, n_p2, n_year] = rowMeans(x_res[index$V, n_p, pts])
          }
        }
      }
      if(output_type == "case+sero" || output_type == "case"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p + n_p0
            output_data$C[n_p2, n_year] = sum(x_res[index$C, n_p, pts])
          }
        }
      }
      if(output_type == "case_alt"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p + n_p0
            output_data$C[, n_p2, n_year] = rowSums(x_res[index$C, n_p, pts])
          }
        }
      }
      if(output_type == "case_alt2"){
        if(n_p0 == 0){
          output_data$day = x_res[1, 1, ]
          output_data$year = x_res[2, 1, ]
        }
        output_data$C[n_p + n_p0, ] = colSums(x_res[index$C, , ])
      }
    }
    x_res = NULL
    gc()
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
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' '
#' @export
#'
parameter_setup <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(), year0 = 1940,
                            vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0, start_SEIRV = list(), mode_time = 0){

  assert_that(mode_start %in% c(0, 1, 2, 3), msg = "mode_start must have value 0, 1, 2 or 3  -  NB 3 should be changed to 1")
  if(mode_start == 3){mode_start = 1} #Temporary fix until mode_start harmonized across all functions/examples
  if(mode_start == 2){assert_that(is.null(start_SEIRV$S) == FALSE, msg = "When mode_start = 2, start_SEIRV data required")}
  assert_that(mode_time %in% c(0:5), msg = "mode_time must be an integer between 0 and 5")
  assert_that(all(FOI_spillover >= 0.0))
  assert_that(all(R0 >= 0.0))
  assert_that(length(pop_data[, 1]) > 1, msg = "Need population data for multiple years")
  assert_that(length(pop_data[1, ]) > 1, msg = "Need population data for multiple age groups")
  n_years = length(pop_data[, 1]) - 1
  N_age = length(pop_data[1, ])
  assert_that(length(vacc_data[, 1]) == n_years + 1, msg = "Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1, ]) == N_age, msg = "No. age groups in population and vaccination data must match")
  assert_that(vaccine_efficacy <= 1.0 && vaccine_efficacy >= 0.0, msg = "Vaccine efficacy must be between 0 and 1")
  assert_that(years_data[1] >= year0, msg = "First data year must be greater than or equal to year0")
  assert_that(max(years_data) + 1 - year0 <= n_years, msg = "Period of years_data must lie within population data")
  assert_that(time_inc %in% c(1, 2.5, 5), msg = "time_inc must have value 1, 2.5 or 5 days")
  pts_year = 365.0/time_inc
  n_t_pts = n_years*pts_year
  n_req = switch(mode_time + 1, 1, n_years, 12, pts_year, n_years*12, n_t_pts)
  assert_that(length(FOI_spillover) == n_req && length(R0) == n_req,
              msg = "Spillover FOI and R0 must be correct length for mode_time")
  inv_365 = 1.0/365.0

  if(mode_time == 0){
    FOI_spillover_t = rep(FOI_spillover, n_t_pts)
    R0_t = rep(R0, n_t_pts)
  } else {
    if(mode_time == 5){
      FOI_spillover_t = FOI_spillover
      R0_t = R0
    } else {
      if(mode_time == 1){
        FOI_spillover_t = R0_t = rep(NA, n_t_pts)
        for(i in 1:n_years){
          FOI_spillover_t[c(1:pts_year) + (i - 1)*pts_year] = rep(FOI_spillover[i], pts_year)
          R0_t[c(1:pts_year) + (i - 1)*pts_year] = rep(R0[i], pts_year)
        }
      } else {
        if(mode_time == 2){date_values = (1 + floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12)}
        if(mode_time == 3){date_values = 1 + (floor(time_inc*c(0:(n_t_pts - 1))) %% pts_year)}
        if(mode_time == 4){date_values = (1 + floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12) + sort(rep(c(1:n_years) - 1,
                                                                                                              pts_year))}
        FOI_spillover_t = FOI_spillover[date_values]
        R0_t = R0[date_values]
      }
    }
  }

  P0 = S_0 = E_0 = I_0 = R_0 = V_0 = rep(0, N_age)
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(N_age, n_years))
  for(i in 1:N_age){ P0[i] = max(1.0, pop_data[1, i]) } #Set all population values to nonzero minimum to avoid NaN values
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i, n_year] = max(1.0, pop_data[n_year + 1, i])*inv_365
      dP2_all[i, n_year] = max(1.0, pop_data[n_year, i])*inv_365
      if(i == 1){
        vacc_rates[i, n_year] = vacc_data[n_year + 1, i]*inv_365
      } else {
        vacc_rates[i, n_year] = max(0.0, vacc_data[n_year + 1, i] - vacc_data[n_year, i - 1])*inv_365
      }
    }
  }

  vacc_initial = vacc_data[1, ]
  if(mode_start == 2){ #Use supplied SEIRV data
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    V_0 = P0*vacc_initial
    if(mode_start == 0){ #No initial immunity
      S_0 = P0*(1.0 - vacc_initial)
    } else { #Stratified herd immunity profile based on notional FOI
      ages = c(1:N_age) - 1
      if(R0_t[1] <= 1.0){
        FOI_estimate = FOI_spillover[1]*365.0
      } else {
        estimation_results = nlm(imm_fraction_function, p =  - 4, R0_t[1], ages, P0/sum(P0))
        FOI_estimate = min(0.1, (FOI_spillover[1]*365.0) + exp(estimation_results$estimate))
      }
      herd_immunity = 1.0 - (exp( - FOI_estimate*(ages + 0.5)))

      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity[i]){
          R_0[i] = P0[i]*(herd_immunity[i] - vacc_initial[i])
          S_0[i] = P0[i]*(1.0 - herd_immunity[i])
        } else {
          S_0[i] = P0[i]*(1.0 - vacc_initial[i])
        }
      }
    }
  }

  return(list(FOI_spillover = FOI_spillover_t, R0 = R0_t, vacc_rate_daily = vacc_rates, N_age = N_age,
              S_0 = S_0, E_0 = E_0, I_0 = I_0, R_0 = R_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
              n_years = n_years, year0 = year0, vaccine_efficacy = vaccine_efficacy, time_inc = time_inc,
              t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious, n_t_pts = n_t_pts))
}
#-------------------------------------------------------------------------------
#' @title imm_fraction_function
#'
#' @description Function to estimate notional FOI for herd immunity based on R0 and population age distribution
#'
#' @details [TBA]
#'
#' @param log_lambda Natural logarithm of force of infection
#' @param R0 Basic reproduction number
#' @param ages List of age values
#' @param pop_fraction Population of each age group as proportion of total
#' '
#' @export
#'
imm_fraction_function <- function(log_lambda =  - 4, R0 = 1.0, ages = c(0:100), pop_fraction = rep(1/101, 101)){
  #TODO  -  Add assert_that functions

  lambda = exp(log_lambda)
  immunity = 1.0 - (exp( - lambda*(ages + 0.5)))
  imm_mean = sum(immunity*pop_fraction)
  imm_mean_target = 1.0 - (1.0/R0)
  dev = abs(imm_mean_target - imm_mean)

  return(dev)
}
#-------------------------------------------------------------------------------
#' @title epi_param_calc
#'
#' @description Calculate FOI_spillover or R0 values from environmental covariates and coefficients
#'
#' @details Takes in environmental covariate values for one or more regions and coefficients of environmental covariates
#'  and calculates epidemiological parameter values via matrix multiplication. Environmental covariates may be constant
#'  or time-varying; constant and time-varying sets of values are supplied as separate input variables, as are coefficients
#'  of constant and time-varying covariates.
#'
#'  The function can accept [TBA]
#'
#' @param coeffs_const Vector of coefficients of time-invariant covariates
#' @param coeffs_var Vector of coefficients of time-varying covariates
#' @param enviro_data_const Data frame of time-invariant environmental covariate values, with region labels in first
#'   column and one row per region
#' @param enviro_data_var List containing data on time-varying environmental covariates: [TBA]
#' '
#' @export
#'
epi_param_calc <- function(coeffs_const = c(), coeffs_var = c(), enviro_data_const = data.frame(), enviro_data_var = NULL){
  #TODO  -  Ensure function works if only variable covariates? (Now works for const only or const + var)
  #TODO  -  Add assertthat checks
  assert_that(is.null(enviro_data_const) == FALSE) #TBC if made capable of using all-variable data
  assert_that(all(c(coeffs_const, coeffs_var) >= 0), msg = "All environmental coefficients must have positive values")
  assert_that(colnames(enviro_data_const)[1] == "region")

  base_output_values = as.vector(as.matrix(enviro_data_const[, c(2:ncol(enviro_data_const))]) %*% as.matrix(coeffs_const))

  if(is.null(enviro_data_var) == FALSE){
    n_pts = dim(enviro_data_var$values)[3]
    var_output_values = colSums(coeffs_var*enviro_data_var$values)
    total_output_values = array(NA, dim = dim(var_output_values))
    for(i in 1:n_pts){ total_output_values[, i] = var_output_values[, i] + base_output_values }
  } else {
    total_output_values = array(base_output_values, dim = c(length(base_output_values), 1))
  }

  return(total_output_values)
}
