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
#' @param seed Random seed (set to NULL if not to be used)
#' '
#' @export
#'
Model_Run2 <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                       year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, output_type = "full", mode_start = 0,
                       start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                       seed = NULL) {

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
                           seed = seed, deterministic = deterministic, preserve_particle_dimension = TRUE)
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
#' @title parameter_setup2
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
parameter_setup2 <- function(FOI_spillover = list(), R0 = list(), vacc_data = list(), pop_data = list(), years_data = c(), year0 = 1940,
                             vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0, start_SEIRV = list(), mode_time = 0){

  #TODO - additional assert_that functions
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
#' @title Generate_Dataset2
#'
#' @description Generate annual serological and/or case/death data
#'
#' @details This function is used to generate annual serological and/or case/death data based on templates;
#' it is normally used by the single_posterior_calc() function.
#'
#' [TBA - Explanation of breakdown of regions to model and how to set lengths of FOI_values and R0_values]
#'
#' @param FOI_values Array of values of force of infection due to spillover from sylvatic reservoir by region + time point
#' @param R0_values Array of values of basic reproduction number for human-human transmission by region and time point
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param sero_template Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#'   and number of samples
#' @param case_template Annual reported case/death data template - data frame with region and year
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_parallel TRUE/FALSE - set model to run in parallel using cluster if TRUE
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param output_frame TRUE/FALSE - indicate whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#' @param seed Optional random seed value to set before running each region for stochastic normalization; set to NULL
#'   to omit; will not work if mode_parallel is not set to FALSE.
#' @param xref_sero Region cross-referencing data generated by template_region_xref() function for serological template;
#'  if set to NULL, calculated within function
#' @param xref_case Region cross-referencing data generated by template_region_xref() function for case template;
#'  if set to NULL, calculated within function
#' '
#' @export
#'
Generate_Dataset2 <- function(FOI_values = c(),R0_values = c(),input_data = list(),sero_template = NULL,case_template = NULL,
                             vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = NULL, mode_time = 0,
                             n_reps = 1,deterministic = FALSE, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                             p_rep_severe = 1.0,p_rep_death = 1.0,mode_parallel = FALSE,cluster = NULL,output_frame = FALSE,
                             seed = NULL, xref_sero=NULL, xref_case=NULL){

  assert_that(input_data_check(input_data),msg = paste("Input data must be in standard format",
                                                       " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template) == FALSE,is.null(case_template) == FALSE),msg = "Need at least one template")
  assert_that(between(vaccine_efficacy,0.0,1.0),msg = "Vaccine efficacy must be between 0-1")
  if(is.null(sero_template) == FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template) == FALSE){ #TODO - add option for monthly data
    assert_that(all(c("region","year") %in% names(case_template)))
    #if("month" %in% names(case_template)){flag_monthly_cases=TRUE}else{flag_monthly_cases=FALSE}
    assert_that(between(p_severe_inf,0.0,1.0),msg = "Severe infection rate must be between 0-1")
    assert_that(between(p_death_severe_inf,0.0,1.0),msg = "Fatality rate of severe infections must be between 0-1")
    assert_that(between(p_rep_severe,0.0,1.0),msg = "Severe infection reporting probability must be between 0-1")
    assert_that(between(p_rep_death,0.0,1.0),msg = "Fatal infection reporting probability must be between 0-1")
  }
  assert_that(is.logical(mode_parallel))
  if(mode_parallel){assert_that(is.null(cluster) == FALSE)}

  #Prune input data based on regions if necessary
  regions = regions_breakdown(c(sero_template$region,case_template$region))
  if(all(regions==input_data$region_labels)==FALSE){
    flag_truncated=TRUE
    input_data = input_data_truncate(input_data,regions)
  } else {
    flag_truncated=FALSE
  }
  n_regions = length(input_data$region_labels)

  #Cross-reference templates with input regions
  if(is.null(sero_template) == FALSE){
    if(is.null(xref_sero) | flag_truncated){xref_sero = template_region_xref(sero_template,input_data$region_labels)}
    sero_line_list = xref_sero$line_list
  } else {
    xref_sero = data.frame(year_data_begin = rep(Inf,n_regions),year_end = rep(-Inf,n_regions))
    sero_line_list = rep(NA,n_regions)
  }
  if(is.null(case_template) == FALSE){
    if(is.null(xref_case) | flag_truncated){xref_case = template_region_xref(case_template,input_data$region_labels)}
    case_line_list = xref_case$line_list
  } else {
    xref_case = data.frame(year_data_begin = rep(Inf,n_regions),year_end = rep(-Inf,n_regions))
    case_line_list = rep(NA,n_regions)
  }
  #TBC
  region_groups=list()
  region_groups[[1]]=which(input_data$region_labels %in% sero_template$region)
  region_groups[[2]]=which(input_data$region_labels %in% case_template$region)
  n_groups=length(region_groups)
  year_data_begin=year_end=rep(NA,n_groups)
  year_data_begin[1]=min(sero_template$year)
  year_end[1]=max(sero_template$year)
  year_data_begin[2]=min(case_template$year)
  year_end[2]=max(case_template$year)

  assert_that(length(dim(FOI_values)) == 2,msg = "FOI_values must be 2-D array")
  assert_that(length(dim(R0_values)) == 2,msg = "R0_values must be 2-D array")
  assert_that(dim(FOI_values)[1] == n_regions,msg = "1st dimension of FOI_values must match number of regions to be modelled")
  assert_that(dim(R0_values)[1] == n_regions,msg = "1st dimension of R0_values must match number of regions to be modelled")
  if(mode_start == 2){assert_that(length(start_SEIRV) == n_regions,
                                  msg = "Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(sero_template)){model_sero_data = NULL} else {
    blank = rep(0,nrow(sero_template))
    model_sero_data = data.frame(samples = blank,positives = blank,sero = blank)
  }
  if(is.null(case_template)){
    model_case_values = model_death_values = NA
  } else {
    model_case_values = model_death_values = rep(0,nrow(case_template))
  }

  #Set up vector of output types to get from model
  #TBC
  output_types=rep(NA,n_groups)
  output_types=c("sero","infs")

  if(mode_parallel){
    #TBA
  }

  for(n_group in 1:n_groups){
    i_regions=region_groups[[n_group]]
    if(mode_parallel){
      #TBA
    } else {
      model_output = Model_Run2(FOI_spillover = array(FOI_values[i_regions,],dim=c(length(i_regions),dim(FOI_values)[2])),
                                R0 = array(R0_values[i_regions,],dim=c(length(i_regions),dim(R0_values)[2])),
                                vacc_data = input_data$vacc_data[i_regions,,],
                                pop_data = input_data$pop_data[i_regions,,],
                                years_data = c(year_data_begin[n_group]:year_end[n_group]),
                                year0 = input_data$years_labels[1], vaccine_efficacy = vaccine_efficacy,
                                time_inc = time_inc, output_type = output_types[n_group],mode_start = mode_start,
                                start_SEIRV = start_SEIRV[[n_group]],mode_time = mode_time,n_particles = n_reps,
                                n_threads = n_reps, deterministic = deterministic, seed = seed)
    }
    t_pts = length(model_output$year)

    for(n_region in i_regions){
      n_region2=match(regions[n_region],regions[region_groups[[n_group]]])
      #Compile case data if needed
      if(is.na(case_line_list[[n_region]][1]) == FALSE){
        case_line_list_region = case_line_list[[n_region]]
        years_case = case_template$year[case_line_list_region]
        n_lines = length(case_line_list_region)

        for(n_rep in 1:n_reps){
          rep_cases = rep_deaths = rep(0,n_lines)
          for(n_line in 1:n_lines){
            pts = c(1:t_pts)[model_output$year == years_case[n_line]]
            infs = sum(model_output$C[n_region2,n_rep,pts])
            if(deterministic){
              severe_infs = floor(infs)*p_severe_inf
              deaths = severe_infs*p_death_severe_inf
              rep_deaths[n_line] = round(deaths*p_rep_death)
              rep_cases[n_line] = rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)
            } else {
              severe_infs = rbinom(1,floor(infs),p_severe_inf)
              deaths = rbinom(1,severe_infs,p_death_severe_inf)
              rep_deaths[n_line] = rbinom(1,deaths,p_rep_death)
              rep_cases[n_line] = rep_deaths[n_line]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
            }
          }
          model_case_values[case_line_list_region] = model_case_values[case_line_list_region]+rep_cases
          model_death_values[case_line_list_region] = model_death_values[case_line_list_region]+rep_deaths
        }
      }

      #Compile seroprevalence data if necessary
      if(is.na(sero_line_list[[n_region]][1]) == FALSE){
        sero_line_list_region = sero_line_list[[n_region]]
        for(n_rep in 1:n_reps){
          for(line in sero_line_list_region){
            #TBC
            i_ages=c((sero_template$age_min[line]+1):sero_template$age_max[line])
            t_pt=match(sero_template$year[line],model_output$year)
            samples=sum(model_output$SEIR_annual[n_region2,i_ages,n_rep,t_pt])
            positives=sum(model_output$R_annual[n_region2,i_ages,n_rep,t_pt])
            model_sero_data$samples[line] = model_sero_data$samples[line]+samples
            model_sero_data$positives[line] = model_sero_data$positives[line] + positives
          }
        }
      }
    }
  }

  if(is.null(sero_template) == FALSE){model_sero_data$sero = model_sero_data$positives/model_sero_data$samples}
  if(is.null(case_template) == FALSE){
    model_case_values = model_case_values/n_reps
    model_death_values = model_death_values/n_reps
  }

  if(output_frame) { #Output complete frames of data
    return(list(model_sero_data = data.frame(region = sero_template$region,year = sero_template$year,
                                             age_min = sero_template$age_min,age_max = sero_template$age_max,
                                             samples = sero_template$samples,
                                             positives = sero_template$samples*model_sero_data$sero,
                                             vc_factor = sero_template$vc_factor),
                model_case_data = data.frame(region = case_template$region,year = case_template$year,
                                             cases = model_case_values,deaths = model_death_values)))
  } else { #Minimal output for MCMC
    return(list(model_sero_values = model_sero_data$sero,model_case_values = model_case_values,
                model_death_values = model_death_values))
  }
}
