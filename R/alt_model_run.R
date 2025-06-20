#TODO:
# -Create new version of Generate_Dataset
# -Run speed tests on desktop and on cluster
#-------------------------------------------------------------------------------
#' @title Model_Run2
#'
#' @description Run SEIRV model for multiple regions as one odin2 run
#'
#' @details TBA
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

  N_age = length(pop_data[1, 1, ]) #Number of age groups
  n_regions = length(pop_data[, 1, 1])
  if(output_type=="full"){
    #model = SEIRV_Model
    step_begin = ((years_data[1] - year0)*(365/time_inc)) #Step at which data starts being saved for final output
    step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1 #Step at which to end
    time_pts = c(step_begin:step_end)
  } else {
    #model = SEIRV_Model2
    i_year_begin=years_data[1] - year0 + 1
    i_year_end=max(years_data) + 1 - year0
    time_pts = (c(i_year_begin:i_year_end)*(365/time_inc))-1
  }
  t_pts_out=length(time_pts)#Number of time points in final output data

  x <- dust_system_create(SEIRV_Model2,
                           pars = parameter_setup2(FOI_spillover, R0,vacc_data, pop_data,
                                                   years_data, year0, vaccine_efficacy, time_inc, mode_start,
                                                   start_SEIRV,mode_time),
                           n_particles = n_particles, n_threads = n_particles, time = 0, dt = 1,
                           deterministic = deterministic, preserve_particle_dimension = TRUE)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = time_pts)
  index = dust_unpack_index(x)

  dim = c(n_regions, N_age, n_particles, t_pts_out)
  if(output_type == "full"){
    dim = c(n_regions, N_age, n_particles, t_pts_out)
    output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                       FOI_total = array(x_res[3, , ]/time_inc, dim = c(n_regions,n_particles,t_pts_out)),
                       S = array(x_res[index$S, , ], dim), E = array(x_res[index$E, , ], dim),
                       I = array(x_res[index$I, , ], dim), R = array(x_res[index$R, , ], dim),
                       V = array(x_res[index$V, , ], dim), C = array(x_res[index$C, , ], dim),
                       C_annual = array(x_res[index$C_annual, , ], dim),
                       R_annual = array(x_res[index$R_annual, , ], dim),
                       SEIR_annual = array(x_res[index$SEIR_annual, , ], dim))

  } else {
    output_data = list(year = years_data)
    if(output_type == "infs"){
      C_annual = array(x_res[index$C_annual, , ], dim)
      output_data$C = array(NA, c(n_regions,n_particles,t_pts_out))
      for(i in 1:n_regions){output_data$C[i,,]=colSums(C_annual[i,,,])}
    } else {
      output_data$R_annual=array(x_res[index$R_annual, , ], dim)
      output_data$SEIR_annual=array(x_res[index$SEIR_annual, , ], dim)
    }
  }
  x_res = NULL
  gc()

  return(output_data)
}
#-------------------------------------------------------------------------------
parameter_setup2 <- function(FOI_spillover = list(), R0 = list(), vacc_data = list(), pop_data = list(), years_data = c(), year0 = 1940,
                             vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0, start_SEIRV = list(), mode_time = 0){

  #TODO - additional assert_that functions
  assert_that(mode_start %in% c(0, 1), msg = "mode_start must have value 0 or 1")
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

  vacc_initial = vacc_data[, 1, ]

  V_0 = P0*vacc_initial
  if(mode_start == 0){ #No initial immunity
    S_0 = P0*(1.0 - vacc_initial)
  } else { #Stratified herd immunity profile based on notional FOI (averaged over first year)
    #TODO
    # R0_year0=mean(R0_t[c(1:pts_year)])
    # FOI_spillover_year0=mean(FOI_spillover_t[c(1:pts_year)])
    # ages = c(1:N_age) - 1
    # if(R0_year0 <= 1.0){
    #   FOI_estimate = FOI_spillover_year0*365.0
    # } else {
    #   estimation_results = nlm(imm_fraction_function, p =  - 4, R0_year0, ages, P0/sum(P0))
    #   FOI_estimate = min(0.1, (FOI_spillover_year0*365.0) + exp(estimation_results$estimate))
    # }
    # herd_immunity = 1.0 - (exp( - FOI_estimate*(ages + 0.5)))
    #
    # for(j in 1:N_age){
    #   if(vacc_initial[i]<herd_immunity[i]){
    #     R_0[i] = P0[i]*(herd_immunity[i] - vacc_initial[i])
    #     S_0[i] = P0[i]*(1.0 - herd_immunity[i])
    #   } else {
    #     S_0[i] = P0[i]*(1.0 - vacc_initial[i])
    #   }
    # }
  }

  return(list(FOI_spillover = FOI_spillover_t, R0 = R0_t, vacc_rate_daily = vacc_rates,
              n_regions = n_regions, N_age = N_age,
              S_0 = S_0, E_0 = E_0, I_0 = I_0, R_0 = R_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
              n_years = n_years, year0 = year0, vaccine_efficacy = vaccine_efficacy, time_inc = time_inc,
              t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious, n_t_pts = n_t_pts))
}
