# R file for general functions in YEP package
#------------------------------------------------
# Global variables
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#------------------------------------------------
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
#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YEP", libpath)
}
#-------------------------------------------------------------------------------
#' @title Model_Run
#'
#' @description Run SEIRV model for single region (Model_Run_Multi_Input can be used to run multiple regions in parallel)
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output:
#'   "full" = SEIRVC + FOI for all steps and ages
#'   "case" = annual total new infections (C) summed across all ages
#'   "sero" = annual SEIRV
#'   "case+sero" = annual SEIRVC, C summed across all ages
#'   "case_alt" = annual total new infections not combined by age
#'   "case_alt2" = total new infections combined by age for all steps
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age,
#'  R0 based only)
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0,
                      vaccine_efficacy = 1.0, time_inc = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
  t_pts_out = step_end-step_begin+1 #Number of time points in final output data

  x <- dust_system_create(SEIRV_Model, pars=parameter_setup(FOI_spillover, R0, vacc_data, pop_data, year0, years_data,
                                                            mode_start,vaccine_efficacy, start_SEIRV, time_inc),
                          n_particles=n_particles,n_threads=n_threads,time=0,dt=1,
                          deterministic=deterministic, preserve_particle_dimension=TRUE)
  index=dust_unpack_index(x)
  dust_system_set_state_initial(x)
  t <- c(step_begin:step_end)
  x_res <- dust_system_simulate(x, t)

  if(output_type == "full"){
    dim = c(N_age, n_particles, t_pts_out)
    output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],FOI_total = x_res[3, , ]/time_inc,
                       S = array(x_res[index$S,,],dim), E = array(x_res[index$E,,],dim), I = array(x_res[index$I,,],dim),
                       R = array(x_res[index$R,,],dim), V = array(x_res[index$V,,],dim), C = array(x_res[index$C,,],dim))
  } else {
    if(output_type == "case_alt2"){
      output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ])
      output_data$C = colSums(x_res[index$C,,])
    }  else {
      n_years = length(years_data)
      output_data = list(year = years_data)
      if(output_type == "case+sero" || output_type == "sero"){
        output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_particles,n_years))
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
          for(n_p in 1:n_particles){
            output_data$C[n_p, n_year] = sum(x_res[index$C, n_p, pts])
          }
        }
      }
      if(output_type == "case_alt"){
        output_data$C = array(0, dim = c(N_age, n_particles, n_years))
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$C[, n_p, n_year] = rowSums(x_res[index$C, n_p, pts])
          }
        }
      }
    }
  }

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
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output:
#'   "full" = SEIRVC + FOI for all steps and ages
#'   "case" = annual total new infections (C) summed across all ages
#'   "sero" = annual SEIRV
#'   "case+sero" = annual SEIRVC, C summed across all ages
#'   "case_alt" = annual total new infections not combined by age
#'   "case_alt2" = total new infections combined by age for all steps
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age,
#'  R0 based only)
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_reps Number of repetitions (used to set number of particles and threads)
#' @param division Number of particles/threads to run in one go (up to 20)
#' '
#' @export
#'
Model_Run_Many_Reps <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(),
                                years_data = c(1940:1941), start_SEIRV = list(), output_type = "full", year0 = 1940,
                                mode_start = 0, vaccine_efficacy = 1.0, time_inc = 1.0, n_reps = 1, division = 10) {

  assert_that(division <= 20, msg = "Number of particles run at once must be 20 or less")
  n_particles0 = min(division, n_reps)
  n_threads = min(division, n_particles0)
  n_divs = ceiling(n_reps/division)
  if(n_divs == 1){
    n_particles_list = n_particles0
  } else {
    n_particles_list = c(rep(n_particles0, n_divs-1), n_reps-(division*(n_divs-1)))
  }

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
  t_pts_out = step_end-step_begin+1 #Number of time points in final output data

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
        output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_reps,n_years))
      }
      if(output_type == "case+sero" || output_type == "case"){
        output_data$C = array(0, dim = c(n_reps, n_years))
      }
      if(output_type == "case_alt"){
        output_data$C = array(0, dim = c(N_age, n_reps, n_years))
      }
    }
  }

  pars=parameter_setup(FOI_spillover, R0, vacc_data, pop_data, year0, years_data,
                       mode_start,vaccine_efficacy, start_SEIRV, time_inc)
  for(div in 1:n_divs){
    n_particles = n_particles_list[div]

    x <- dust_system_create(SEIRV_Model, pars=pars,
                            n_particles=n_particles,n_threads=n_threads,time=0,dt=1,
                            deterministic=FALSE, preserve_particle_dimension=TRUE)
    dust_system_set_state_initial(x)
    t <- c(step_begin:step_end)
    x_res <- dust_system_simulate(x, t)

    if(div == 1){
      n_p0 = 0
      index=dust_unpack_index(x)
    }else{n_p0 = sum(n_particles_list[c(1:(div-1))])}

    if(output_type == "full"){
      n_p_values = c(1:n_particles)+n_p0
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
            n_p2 = n_p+n_p0
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
            n_p2 = n_p+n_p0
            output_data$C[n_p2, n_year] = sum(x_res[index$C, n_p, pts])
          }
        }
      }
      if(output_type == "case_alt"){
        for(n_year in 1:n_years){
          pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2 = n_p+n_p0
            output_data$C[, n_p2, n_year] = rowSums(x_res[index$C, n_p, pts])
          }
        }
      }
      if(output_type == "case_alt2"){
        if(n_p0 == 0){
          output_data$day = x_res[1, 1, ]
          output_data$year = x_res[2, 1, ]
        }
        output_data$C[n_p+n_p0,] = colSums(x_res[index$C,,])
        # for(pt in 1:t_pts_out){
        #   for(n_p in 1:n_particles){
        #     n_p2 = n_p+n_p0
        #     output_data$C[n_p2, pt] = sum(x_res[index$C, n_p, pt])
        #   }
        # }
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
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param year0 First year in population/vaccination data
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age,
#'  R0 based only)
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' '
#' @export
#'
parameter_setup <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), year0 = 1940,
                            years_data = c(1941:1942), mode_start = 0, vaccine_efficacy = 1.0, start_SEIRV = list(),
                            time_inc = 1.0){

  assert_that(FOI_spillover>=0.0)
  assert_that(R0>=0.0)
  assert_that(length(pop_data[, 1])>1, msg = "Need population data for multiple years")
  assert_that(length(pop_data[1, ])>1, msg = "Need population data for multiple age groups")
  n_years = length(pop_data[, 1])-1
  N_age = length(pop_data[1, ])
  assert_that(length(vacc_data[, 1]) == n_years+1, msg = "Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1, ]) == N_age, msg = "No. age groups in population and vaccination data must match")
  assert_that(mode_start %in% c(0, 1, 2, 3), msg = "mode_start must have value 0, 1, 2 or 3")
  assert_that(vaccine_efficacy <= 1.0 && vaccine_efficacy >= 0.0, msg = "Vaccine efficacy must be between 0 and 1")
  if(mode_start == 2){
    assert_that(is.null(start_SEIRV$S) == FALSE, msg = "When mode_start = 2, start_SEIRV data is required")
  }
  assert_that(years_data[1] >= year0, msg = "First data year must be greater than or equal to year0")
  assert_that(max(years_data)+1-year0 <= n_years, msg = "Period of years_data must lie within population data")
  vacc_initial = vacc_data[1, ]
  assert_that(time_inc %in% c(1, 2.5, 5),
              msg = "time_inc must have value 1, 2.5 or 5 days (must have integer no. points/year)")
  inv_365 = 1.0/365.0

  P0 = S_0 = E_0 = I_0 = R_0 = V_0 = rep(0, N_age)
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(N_age, n_years))
  for(i in 1:N_age){
    P0[i] = max(1.0, pop_data[1, i]) #Set all population values to nonzero minimum to avoid NaN values
  }
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i, n_year] = max(1.0, pop_data[n_year+1, i])*inv_365
      dP2_all[i, n_year] = max(1.0, pop_data[n_year, i])*inv_365
      if(i == 1){
        vacc_rates[i, n_year] = vacc_data[n_year+1, i]*inv_365
      } else {
        vacc_rates[i, n_year] = max(0.0, vacc_data[n_year+1, i]-vacc_data[n_year, i-1])*inv_365
      }
    }
  }

  #-----------------------------------------------------------------------------
  if(mode_start == 2){
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    V_0 = P0*vacc_initial
    #-----------------------------------------------------------------------------
    if(mode_start == 0){
      S_0 = P0*(1.0-vacc_initial)
    }
    #-----------------------------------------------------------------------------
    if(mode_start == 1){ #Herd immunity, uniform by age
      if(R0>1.0){
        herd_immunity = 1.0-(1.0/R0)
      } else {
        herd_immunity = 0.0
      }
      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity){
          R_0[i] = P0[i]*(herd_immunity-vacc_initial[i])
          S_0[i] = P0[i]*(1.0-herd_immunity)
        } else {
          S_0[i] = P0[i]*(1.0-vacc_initial[i])
        }
      }
    }
    #-----------------------------------------------------------------------------
    if(mode_start == 3){ #New herd immunity calculation to give age-stratified immunity profile based on notional FOI
      ages = c(1:N_age)-1
      if(R0 <= 1.0){
        FOI_estimate = FOI_spillover*365.0
      } else {
        estimation_results = nlm(imm_fraction_function, p = -4, R0, ages, P0/sum(P0))
        FOI_estimate = min(0.1, (FOI_spillover*365.0)+exp(estimation_results$estimate))
      }
      herd_immunity = 1.0-(exp(-FOI_estimate*(ages+0.5)))

      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity[i]){
          R_0[i] = P0[i]*(herd_immunity[i]-vacc_initial[i])
          S_0[i] = P0[i]*(1.0-herd_immunity[i])
        } else {
          S_0[i] = P0[i]*(1.0-vacc_initial[i])
        }
      }
    }
  }

  return(list(FOI_spillover = FOI_spillover, R0 = R0, vacc_rate_daily = vacc_rates, N_age = N_age,
              S_0 = S_0, E_0 = E_0, I_0 = I_0, R_0 = R_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
              n_years = n_years, year0 = year0, vaccine_efficacy = vaccine_efficacy, time_inc = time_inc,
              t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious))
}
#-------------------------------------------------------------------------------
#' @title create_param_labels
#'
#' @description Apply names to the parameters in a set used for data matching and parameter fitting
#'
#' @details Takes in environmental covariate data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param enviro_data Environmental data frame, containing only relevant environmental covariate values for regions of
#'                    interest
#' @param extra_estimated_params Vector of strings listing variable parameters besides ones determining FOI/R0 (may include
#' vaccine efficacy and/or infection/death reporting probabilities and/or Brazil FOI adjustment factor)
#'
#' @export
#'
create_param_labels <- function(enviro_data = NULL, extra_estimated_params = c("vacc_eff")){

  assert_that(is.data.frame(enviro_data), msg = "Environmental covariate values must be in data frame")

  n_extra = length(extra_estimated_params)
  env_vars = colnames(enviro_data)[c(2:ncol(enviro_data))]
  n_env_vars = length(env_vars)
  n_params = (2*n_env_vars)+n_extra
  param_names = rep("", n_params)
  for(i in 1:n_env_vars){
    param_names[i] = paste("FOI_", env_vars[i], sep = "")
    param_names[i+n_env_vars] = paste("R0_", env_vars[i], sep = "")
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params] = extra_estimated_params}

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
param_calc_enviro <- function(enviro_coeffs = c(), enviro_covar_values = c()){

  assert_that(all(enviro_coeffs >= 0), msg = "All environmental coefficients must have positive values")
  n_env_vars = length(enviro_covar_values)
  assert_that(length(enviro_coeffs)==2*n_env_vars, msg = "Wrong number of environmental coefficients")

  output = list()
  output$FOI = sum(enviro_coeffs[c(1:n_env_vars)]*enviro_covar_values)
  output$R0 = sum(enviro_coeffs[c(1:n_env_vars)+n_env_vars]*enviro_covar_values)

  return(output)
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
imm_fraction_function <- function(log_lambda = -4, R0 = 1.0, ages = c(0:100), pop_fraction = rep(1/101, 101)){
  #TODO - Add assert_that functions

  lambda = exp(log_lambda)
  immunity = 1.0-(exp(-lambda*(ages+0.5)))
  imm_mean = sum(immunity*pop_fraction)

  imm_mean_target = 1.0-(1.0/R0)
  dev = abs(imm_mean_target-imm_mean)

  return(dev)
}
#-------------------------------------------------------------------------------
#Function not currently in use - potential way to run for multiple regions at once, possibly not improvement, test in odin2
#' #' @title Model_Run_Multi_Input
#' #'
#' #' @description Run SEIRV model for multiple parameter sets
#' #'
#' #' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' #' for multiple parameter sets (spillover FOI, R0, population and vaccination data, starting data values where relevant)
#' #' over a specified time period for a number of particles/threads and outputs time-dependent SEIRV values, infection
#' #' numbers and/or total force of infection values. Parameter sets may represent multiple regions and/or the same region
#' #' with different inputs.
#' #'
#' #' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#' #' @param R0 Vector of values of basic reproduction number for urban spread of infection
#' #' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by parameter set,
#'                     age group and year
#' #' @param pop_data Population by parameter set and age group by year
#' #' @param years_data Incremental vector of years denoting years for which to save data
#' #' @param start_SEIRV SEIRV data from end of a previous run to use as input, for each parameter set
#' #' @param output_type Type of data to output:
#' #'   "full" = SEIRVC + FOI for all steps and ages
#' #'   "case" = annual total new infections (C) summed across all ages
#' #'   "sero" = annual SEIRV
#' #'   "case+sero" = annual SEIRVC, C summed across all ages
#' #'   "case_alt" = annual total new infections not combined by age
#' #'   "case_alt2" = total new infections combined by age for all steps
#' #' @param year0 First year in population/vaccination data
#' #' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#' #'  If mode_start = 0, only vaccinated individuals
#' #'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0
#'      based only)
#' #'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' #'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' #' @param vaccine_efficacy Proportional vaccine efficacy
#' #' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' #' @param n_particles number of particles to use (up to 20)
#' #' @param n_threads number of threads to use
#' #' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' #' '
#' #' @export
#' #'
#' Model_Run_Multi_Input <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(),
#'                                   years_data = c(1940:1941), start_SEIRV = list(), output_type = "full", year0 = 1940,
#'                                   mode_start = 0, vaccine_efficacy = 1.0, time_inc = 1.0, n_particles = 1, n_threads=1,
#'                                   deterministic = FALSE) {
#'
#'   #TODO Add assert_that functions
#'   assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")
#'
#'   n_nv = 3 #Number of non-vector outputs
#'   N_age = length(pop_data[1, 1, ]) #Number of age groups
#'   n_data_pts = (6*N_age)+n_nv #Number of data values per time point in output
#'   step_begin = ((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
#'   step_end = ((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
#'   t_pts_out = step_end-step_begin+1 #Number of time points in final output data
#'   #dim = c(N_age, n_particles, t_pts_out)
#'
#'   n_param_sets = length(FOI_spillover)
#'   pars = list()
#'   for(i in 1:n_param_sets){
#'     pars[[i]] = parameter_setup(FOI_spillover[i], R0[i], vacc_data[i, , ], pop_data[i, , ], year0, years_data,
#'                                 mode_start, vaccine_efficacy, start_SEIRV[[i]], time_inc)
#'   }
#'
#'   x <- SEIRV_Model$new(pars, time = 1, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic,
#'                        pars_multi = TRUE)
#'
#'   x_res <- array(NA, dim = c(n_data_pts, n_particles*n_param_sets, t_pts_out))
#'   for(step in step_begin:step_end){
#'     x_res[, , step-step_begin+1] <- x$run(step)
#'   }
#'   if(step_begin == 0){x_res[2, , 1] = rep(year0, n_particles*n_param_sets)}
#'
#'   output = list()
#'   for(i in 1:n_param_sets){
#'     n_p2 = c(1:n_particles)+(n_particles*(i-1))
#'     if(output_type == "full"){
#'       dim = c(N_age, n_particles, t_pts_out)
#'       output[[i]] = list(day = x_res[1, n_p2[1], ], year = x_res[2, n_p2[1], ])
#'       output[[i]]$FOI_total = array(x_res[3, n_p2, ]/time_inc, dim = c(n_particles, t_pts_out))
#'       output[[i]]$S = array(x_res[index$S, n_p2, ], dim)
#'       output[[i]]$E = array(x_res[index$E, n_p2, ], dim)
#'       output[[i]]$I = array(x_res[index$I, n_p2, ], dim)
#'       output[[i]]$R = array(x_res[index$R, n_p2, ], dim)
#'       output[[i]]$V = array(x_res[index$V, n_p2, ], dim)
#'       output[[i]]$C = array(x_res[index$C, n_p2, ], dim)
#'     } else {
#'       if(output_type == "case_alt2"){
#'         output[[i]] = list(day = x_res[1, 1, ], year = x_res[2, 1, ])
#'         output[[i]]$C = array(0, dim = c(n_particles, t_pts_out))
#'         for(pt in 1:t_pts_out){
#'           for(n_p in 1:n_particles){
#'             output[[i]]$C[n_p, pt] = sum(x_res[index$C, n_p, pt])
#'           }
#'         }
#'       } else {
#'         n_years = length(years_data)
#'         dim = c(N_age, n_particles, n_years)
#'         output[[i]] = list(year = years_data)
#'         if(output_type == "case+sero" || output_type == "sero"){
#'           output[[i]]$V = output[[i]]$R = output[[i]]$I = output[[i]]$E = output[[i]]$S = array(0, dim)
#'           for(n_year in 1:n_years){
#'             pts = c(1:t_pts_out)[x_res[2, n_p2[1], ] == years_data[n_year]]
#'             for(n_p in 1:n_particles){
#'               output[[i]]$S[, n_p, n_year] = rowMeans(x_res[index$S, n_p2[n_p], pts])
#'               output[[i]]$E[, n_p, n_year] = rowMeans(x_res[index$E, n_p2[n_p], pts])
#'               output[[i]]$I[, n_p, n_year] = rowMeans(x_res[index$I, n_p2[n_p], pts])
#'               output[[i]]$R[, n_p, n_year] = rowMeans(x_res[index$R, n_p2[n_p], pts])
#'               output[[i]]$V[, n_p, n_year] = rowMeans(x_res[index$V, n_p2[n_p], pts])
#'             }
#'           }
#'         }
#'         if(output_type == "case+sero" || output_type == "case"){
#'           output[[i]]$C = array(0, dim = c(n_particles, n_years))
#'           for(n_year in 1:n_years){
#'             pts = c(1:t_pts_out)[x_res[2, n_p2[1], ] == years_data[n_year]]
#'             for(n_p in 1:n_particles){
#'               output[[i]]$C[n_p, n_year] = sum(x_res[index$C, n_p2[n_p], pts])
#'             }
#'           }
#'         }
#'         if(output_type == "case_alt"){
#'           output[[i]]$C = array(0, dim = c(N_age, n_particles, n_years))
#'           for(n_year in 1:n_years){
#'             pts = c(1:t_pts_out)[x_res[2, 1, ] == years_data[n_year]]
#'             for(n_p in 1:n_particles){
#'               output[[i]]$C[, n_p, n_year] = rowSums(x_res[index$C, n_p, pts])
#'             }
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   return(output)
#' }
