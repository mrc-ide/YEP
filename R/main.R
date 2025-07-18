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
#' @importFrom dplyr between
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
# TODO - Add provision for monthly case data
# TODO - Update documentation
# TODO - Add provision for additional output data formats like original (in progress)
#' @title Model_Run
#'
#' @description Run SEIRV model for one or more regions using odin2 models
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model for one or
#' more regions over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values. Run up to 20 particles (Model_Run_Many_Reps
#' can be used to run larger numbers of particles).
#'
#' @param FOI_spillover Matrix of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Matrix of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by region. age group and year
#' @param pop_data Population by region, age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_out Type of data to output: \cr
#'   If mode_out = 1, SEIRVC + FOI for all regions, steps and ages at every time point in years_data \cr
#'   If mode_out = 2 annual totals of SEIR, R and V by region and age for calculating seroprevalence
#'   If mode_out = 3, annual total new infections (C_annual) by region combined across age groups \cr
#'   If mode_out = 4, annual total new infections (C_annual) not combined by age \cr
#'   If mode_out = 5, total new infections (C) combined by age at every time point in years_data \cr
#'   If mode_out = 6 [TBA - combine outputs for 2 and 3, for case+sero data]
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) (TBD) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' @param seed Random seed (set to NULL if not to be used)
#' '
#' @export
#'
Model_Run <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_out = 1, mode_start = 0,
                      start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                      seed = NULL) {

  #TODO Add assert_that functions? (NB  -  Some checks carried out in parameter_setup)assert_that(mode_out %in% c(1:5))
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, 1, ]) #Number of age groups
  n_regions = length(pop_data[, 1, 1])
  model = switch(mode_out,
                 SEIRV_Model_mr01_basic,SEIRV_Model_mr02_sero,SEIRV_Model_mr03_infs,
                 SEIRV_Model_mr03_infs,SEIRV_Model_mr01_basic)
  if(mode_out %in% c(2:4)){
    i_year_begin=years_data[1] - year0 + 1
    i_year_end=max(years_data) + 1 - year0
    time_pts = (c(i_year_begin:i_year_end)*(365/time_inc))-1
  } else {
    step_begin = ((years_data[1] - year0)*(365/time_inc))
    step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1
    time_pts = c(step_begin:step_end)
  }
  t_pts_out=length(time_pts) #Number of time points in final output data

  x <- dust_system_create(model,
                          pars = parameter_setup(FOI_spillover, R0,vacc_data, pop_data,
                                                 years_data, year0, vaccine_efficacy,
                                                 time_inc, mode_start, start_SEIRV, mode_time),
                          n_particles = n_particles, n_threads = n_particles, time = 0,
                          dt = 1, seed = seed, deterministic = deterministic,
                          preserve_particle_dimension = TRUE)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = time_pts)
  index = dust_unpack_index(x)

  dim1 = c(n_regions, N_age, n_particles, t_pts_out)
  if(mode_out %in% c(2,3)){
    output_data = list(year = years_data)
    if(mode_out == 2){
      output_data$R_annual=array(x_res[index$R_annual, , ], dim1)
      output_data$SEIR_annual=array(x_res[index$SEIR_annual, , ], dim1)
      output_data$V_annual=array(x_res[index$V_annual, , ], dim1)
    } else {
      C_all = array(x_res[index$C_annual, , ], dim1)
      output_data$C_annual = array(0, dim = c(n_regions,n_particles,t_pts_out))
      for(i in 1:n_regions){output_data$C_annual[i,,]=colSums(C_all[i,,,])}
    }
  } else {
    if(mode_out == 1){
      output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                         FOI_total = array(x_res[3, , ]/time_inc, dim = c(n_regions,n_particles,t_pts_out)),
                         S = array(x_res[index$S, , ], dim1), E = array(x_res[index$E, , ], dim1),
                         I = array(x_res[index$I, , ], dim1), R = array(x_res[index$R, , ], dim1),
                         V = array(x_res[index$V, , ], dim1), C = array(x_res[index$C, , ], dim1))
    } else {
      if(mode_out == 4){
        output_data = list(year = years_data,
                           C_annual = array(x_res[index$C_annual, , ], dim1))
      }
      if(mode_out == 5){
        output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                           C = array(0, dim = c(n_regions, n_particles, t_pts_out)))
        C_all=array(x_res[index$C, , ], dim1)
        for(i in 1:n_regions){output_data$C[i,,]=colSums(C_all[i,,,])}
      }
    }
  }
  x_res = NULL
  gc()

  return(output_data)
}
#-------------------------------------------------------------------------------
# TODO - Adapt to new model version
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
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_out Type of data to output: \cr
#'   "full" = SEIRVC + FOI for all steps and ages \cr
#'   "infs" = annual total new infections (C) summed across all ages \cr
#'   "sero" = annual SEIRV \cr
#'   "infs_sero" = annual SEIRV, C summed across all ages \cr
#'   "infs_alt" = annual total new infections not combined by age \cr
#'   "infs_alt2" = total new infections combined by age for all steps
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' @param n_reps Number of repetitions (used to set number of particles and threads)
#' @param division Number of particles/threads to run in one go (up to 20)
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Many_Reps <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                                year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_out = "full", mode_start = 0,
                                start_SEIRV = list(), mode_time = 0, n_reps = 1, division = 10, deterministic = FALSE) {

  assert_that(division <= 20, msg = "Number of particles run at once must be 20 or less")
  assert_that(mode_out %in% c("full","infs","sero","infs_sero","infs_alt","infs_alt2"))
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

  if(mode_out == "full"){
    dim = c(N_age, n_reps, t_pts_out)
    output_data = list(day = rep(NA, t_pts_out), year = rep(NA, t_pts_out), FOI_total = array(NA, c(n_reps, t_pts_out)),
                       S = array(NA, dim), E = array(NA, dim), I = array(NA, dim),
                       R = array(NA, dim), V = array(NA, dim), C = array(NA, dim))
  } else {
    if(mode_out == "infs_alt2"){
      output_data = list(day = rep(NA, t_pts_out), year = rep(NA, t_pts_out))
      output_data$C = array(0, dim = c(n_reps, t_pts_out))
    } else {
      n_years = length(years_data)
      output_data = list(year = years_data)
      if(mode_out == "infs_sero" || mode_out == "sero"){
        output_data$V = output_data$R = output_data$I = output_data$E = output_data$S = array(0, dim = c(N_age, n_reps, n_years))
      }
      if(mode_out == "infs_sero" || mode_out == "infs"){ output_data$C = array(0, dim = c(n_reps, n_years)) }
      if(mode_out == "infs_alt"){ output_data$C = array(0, dim = c(N_age, n_reps, n_years)) }
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

    if(mode_out == "full"){
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
      if(mode_out == "infs_sero" || mode_out == "sero"){
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
      if(mode_out == "infs_sero" || mode_out == "infs"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p + n_p0
            output_data$C[n_p2, n_year] = sum(x_res[index$C, n_p, pts])
          }
        }
      }
      if(mode_out == "infs_alt"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p + n_p0
            output_data$C[, n_p2, n_year] = rowSums(x_res[index$C, n_p, pts])
          }
        }
      }
      if(mode_out == "infs_alt2"){
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
#' @title parameter_setup
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin SEIRV model.
#'
#' @param FOI_spillover Matrix of values of force of infection due to spillover from sylvatic reservoir TBA
#' @param R0 Matrix of values of basic reproduction number for urban spread of infection TBA
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by region, age group and year
#' @param pop_data Population by region, age group and year
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
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' '
#' @export
#'
parameter_setup <- function(FOI_spillover = list(), R0 = list(), vacc_data = list(), pop_data = list(),
                            years_data = c(), year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0,
                            start_SEIRV = list(), mode_time = 0){

  #TODO - additional assert_that functions?
  assert_that(mode_start %in% c(0, 1, 2), msg = "mode_start must have value 0, 1 or 2")
  assert_that(mode_time %in% c(0:5), msg = "mode_time must be an integer between 0 and 5")
  assert_that(all(FOI_spillover >= 0.0))
  assert_that(all(R0 >= 0.0))
  assert_that(length(pop_data[1, , 1]) > 1, msg = "Need population data for multiple years")
  assert_that(length(pop_data[1, 1, ]) > 1, msg = "Need population data for multiple age groups")
  n_regions = length(pop_data[, 1, 1])
  n_years = length(pop_data[1, , 1]) - 1
  N_age = length(pop_data[1, 1, ])
  assert_that(length(vacc_data[1, , 1]) == n_years + 1,
              msg = "Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[, 1, 1])==n_regions, msg = "No. regions in population and vaccination data must match")
  assert_that(length(vacc_data[1, 1, ]) == N_age, msg = "No. age groups in population and vaccination data must match")
  assert_that(between(vaccine_efficacy, 0.0, 1.0), msg = "Vaccine efficacy must be between 0-1")
  assert_that(years_data[1] >= year0, msg = "First data year must be greater than or equal to year0")
  assert_that(max(years_data) + 1 - year0 <= n_years, msg = "Period of years_data must lie within population data")
  assert_that(time_inc %in% c(1, 2.5, 5), msg = "time_inc must have value 1, 2.5 or 5 days")
  pts_year = 365.0/time_inc
  n_t_pts = n_years*pts_year
  n_req = switch(mode_time + 1, 1, n_years, 12, pts_year, n_years*12, n_t_pts)
  assert_that(dim(FOI_spillover)[2] == n_req && dim(R0)[2] == n_req,
              msg = "Spillover FOI and R0 must be correct length for mode_time")
  inv_365 = 1.0/365.0

  date_values = switch(mode_time + 1,
                       rep(1, n_t_pts),
                       sort(rep(c(1:n_years), pts_year)),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12),
                       1 + (floor(time_inc*c(0:(n_t_pts - 1))) %% pts_year),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12) + (12*sort(rep(c(1:n_years) - 1,
                                                                                                pts_year))),
                       c(1:n_t_pts))

  FOI_spillover_t = array(FOI_spillover[, date_values],dim=c(n_regions,n_t_pts))
  R0_t = array(R0[, date_values],dim=c(n_regions,n_t_pts))

  P0 = S_0 = E_0 = I_0 = R_0 = V_0 = array(0, dim = c(n_regions, N_age))
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(n_regions, N_age, n_years))
  for(i in 1:n_regions){for(j in 1:N_age){ P0[i,j] = max(1.0, pop_data[i,1, j]) }}
  for(n_year in 1:n_years){
    for(i in 1:n_regions){
      for(j in 1:N_age){
        dP1_all[i,j,n_year] = max(1.0, pop_data[i, n_year + 1, j])*inv_365
        dP2_all[i,j,n_year] = max(1.0, pop_data[i, n_year, j])*inv_365
        if(j == 1){
          vacc_rates[i,j,n_year] = vacc_data[i, n_year + 1, j]*inv_365
        } else {
          vacc_rates[i,j,n_year] = max(0.0, vacc_data[i, n_year + 1, j] - vacc_data[i, n_year, j - 1])*inv_365
        }
      }
    }
  }

  if(mode_start == 2){ #Use supplied SEIRV data
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    vacc_initial = array(vacc_data[, 1, ],dim=c(n_regions,N_age))
    V_0 = P0*vacc_initial
    if(mode_start == 0){ #No initial immunity
      S_0 = P0*(1.0 - vacc_initial)
    } else { #Stratified herd immunity profile based on notional FOI (averaged over first year)
      ages = c(1:N_age) - 1
      for(i in 1:n_regions){
        R0_year0=mean(R0_t[i,c(1:pts_year)])
        FOI_spillover_year0=mean(FOI_spillover_t[i,c(1:pts_year)])
        if(R0_year0 <= 1.0){
          FOI_estimate = FOI_spillover_year0*365.0
        } else {
          estimation_results = nlm(imm_fraction_function, p =  - 4, R0_year0, ages, P0[i,]/sum(P0[i,]))
          FOI_estimate = min(0.1, (FOI_spillover_year0*365.0) + exp(estimation_results$estimate))
        }
        herd_immunity = 1.0 - (exp( - FOI_estimate*(ages + 0.5)))

        for(j in 1:N_age){
          if(vacc_initial[i,j]<herd_immunity[j]){
            R_0[i,j] = P0[i,j]*(herd_immunity[j] - vacc_initial[i,j])
            S_0[i,j] = P0[i,j]*(1.0 - herd_immunity[j])
          } else {
            S_0[i,j] = P0[i,j]*(1.0 - vacc_initial[i,j])
          }
        }
      }
    }
  }

  return(list(FOI_spillover = FOI_spillover_t, R0 = R0_t, vacc_rate_daily = vacc_rates,
              n_regions = n_regions, N_age = N_age,
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
epi_param_calc <- function(coeffs_const = c(0), coeffs_var = c(0), enviro_data_const = data.frame(), enviro_data_var = NULL){
  #TODO  -  Ensure function works if only variable covariates? (Now works for const only or const + var)
  #TODO  -  Add assertthat checks
  assert_that(is.null(enviro_data_const) == FALSE) #TBC if made capable of using all-variable data
  assert_that(all(c(coeffs_const, coeffs_var) >= 0), msg = "All environmental coefficients must have positive values")
  assert_that(colnames(enviro_data_const)[1] == "region", msg = "Constant environmental data must contain regions")
  assert_that(is.numeric(coeffs_const) && is.numeric(coeffs_var),msg="Coefficients must be in numerical vectors")

  base_output_values = as.vector(as.matrix(enviro_data_const[, c(2:ncol(enviro_data_const))]) %*% as.matrix(coeffs_const))

  if(is.null(enviro_data_var) == FALSE){
    #assert_that(all(enviro_data_const$region))
    n_pts = dim(enviro_data_var$values)[3]
    var_output_values = colSums(coeffs_var*enviro_data_var$values)
    total_output_values = array(NA, dim = dim(var_output_values))
    for(i in 1:n_pts){ total_output_values[, i] = var_output_values[, i] + base_output_values }
  } else {
    total_output_values = array(base_output_values, dim = c(length(base_output_values), 1))
  }

  return(total_output_values)
}

