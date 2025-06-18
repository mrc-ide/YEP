#' @title Model_Run2
#'
#' @description Run SEIRV model for single region (alternate version)
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
#'   "infs" = annual total new infections (C) summed across all ages \cr
#'   "sero" = TBA
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
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run2 <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, output_type = "full", mode_start = 0,
                      start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB  -  Some checks carried out in parameter_setup)
  assert_that(output_type %in% c("full","infs","sero"))
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  if(output_type=="full"){
    step_begin = ((years_data[1] - year0)*(365/time_inc)) #Step at which data starts being saved for final output
    step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1 #Step at which to end
    time_pts = c(step_begin:step_end)
  } else {
    i_year_begin=years_data[1] - year0 + 1
    i_year_end=max(years_data) + 1 - year0
    time_pts = (c(i_year_begin:i_year_end)*(365/time_inc))-1
  }
  t_pts_out=length(time_pts)#Number of time points in final output data

  x <- dust_system_create(SEIRV_Model2, pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                                                              vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time),
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = time_pts)

  if(output_type == "full"){
    dim = c(N_age, n_particles, t_pts_out)
    output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ], FOI_total = x_res[3, , ]/time_inc,
                       S = array(x_res[index$S, , ], dim), E = array(x_res[index$E, , ], dim), I = array(x_res[index$I, , ], dim),
                       R = array(x_res[index$R, , ], dim), V = array(x_res[index$V, , ], dim), C = array(x_res[index$C, , ], dim),
                       C_annual = array(x_res[index$C_annual, , ], dim))

  } else {
    output_data = list(year = years_data)
    if(output_type == "infs"){
      output_data$C = array(colSums(x_res[index$C_annual, , ]), dim=c(n_particles,t_pts_out))
    } else {
      dim = c(N_age, n_particles, t_pts_out)
      output_data$R_annual=array(x_res[index$R_annual, , ], dim)
      output_data$SEIR_annual=array(x_res[index$SEIR_annual, , ], dim)
    }
  }
  x_res = NULL
  gc()

  return(output_data)
}
