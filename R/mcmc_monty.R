# File for functions used to estimate parameters via Monty
# TODO:
# -Add check function(s) and/or assert_that checks in existing functions
# -Make separate function for new region indexing?
#-------------------------------------------------------------------------------
#' @title pars_var_setup
#'
#' @description Set up variable parameters for fitting
#'
#' @details Construct data frame of parameters to be estimated via MCMC, containing
#'  names, maximum and minimum permitted values, and mean and standard deviation
#'  for prior calculations
#'
#' @param n_env_vars Number of environmental covariates
#' @param vars_extra_names Names of additional parameters, out of:\cr
#'  "vaccine_efficacy","p_severe_inf","p_death_severe_inf","p_rep_severe",\cr
#'  "p_rep_death","m_FOI_BRA", "a_T0", "a_Tm", "a_c", "mu_T0", "mu_Tm", \cr
#'  "mu_c","PDR_T0", "PDR_Tm","PDR_c","log_a_ptf", "overdisp"\cr
#'  (Definitions TBA)
#'
#' @export
#'
pars_var_setup <- function(n_env_vars = 5,vars_extra_names = c("p_rep_severe","p_rep_death")){

  assert_that(all(vars_extra_names %in% extra_param_names))
  n_extra = length(vars_extra_names)
  n_rows = n_extra + (2*n_env_vars)
  pars_var = data.frame(name = c(vars_extra_names,paste0("log_FOI_coeffs[",c(1:n_env_vars),"]"),
                                 paste0("log_R0_coeffs[",c(1:n_env_vars),"]")),
                        max = c(rep(1,n_extra),rep(-10,n_env_vars),rep(1,n_env_vars)),
                        min = c(rep(0.05,n_extra),rep(-20,n_env_vars),rep(-10,n_env_vars)),
                        mean = c(rep(1,n_extra),rep(-15,n_env_vars),rep(-3,n_env_vars)),
                        sd = rep(1000,n_rows))

  return(pars_var)
}
#-------------------------------------------------------------------------------
#' @title pars_fixed_setup
#'
#' @description Set up fixed parameters for fitting
#'
#' @details TBA
#'
#' @param sero_template Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param case_template Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by region. age group and year
#' @param pop_data Population by region, age group and year
#' @param years_data Vector of years denoting years for which data needed (TBC?)
#' @param year0 First year in population/vaccination data
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) (TBD) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param fixed_extra List containing additional fixed parameters (from list available\cr
#'    as variable parameters)
#' @param ref_BRA List of region numbers for which Brazil FOI multiplier to be applied
#' @param i_ttf Number of environmental covariate for which temperature transformation to be applied
#' @param i_ptf Number of environmental covariate for which precipitation transformation to be applied
#'
#' @export
#'
pars_fixed_setup <- function(sero_template = list(),case_template = list(), vacc_data = list(),
                             pop_data = list(), years_data = c(), year0 = 1940,  time_inc = 1.0,
                             mode_start = 0, start_SEIRV = NULL, fixed_extra = list(), ref_BRA = NULL,
                             i_ttf = NULL, i_ptf = NULL){

  assert_that(length(pop_data[1, , 1]) > 1, msg = "Need population data for multiple years")
  assert_that(length(pop_data[1, 1, ]) > 1, msg = "Need population data for multiple age groups")
  n_regions = length(pop_data[, 1, 1])
  n_years = length(pop_data[1, , 1]) - 1
  N_age = length(pop_data[1, 1, ])
  assert_that(length(vacc_data[1, , 1]) ==  n_years + 1,
              msg = "Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[, 1, 1]) == n_regions, msg = "No. regions in population and vaccination data must match")
  assert_that(length(vacc_data[1, 1, ]) ==  N_age, msg = "No. age groups in population and vaccination data must match")
  assert_that(years_data[1] >=  year0, msg = "First data year must be greater than or equal to year0")
  #Need 2 years of leeway - TBA
  assert_that(max(years_data) + 2 - year0 <=  n_years, msg = "Period of years_data must lie within population data with 2 years at end")
  assert_that(time_inc %in% c(1, 2.5, 5), msg = "time_inc must have value 1, 2.5 or 5 days")
  assert_that(mode_start %in% c(0, 1, 2), msg = "mode_start must have value 0, 1 or 2")
  if(mode_start == 2){assert_that(is.null(start_SEIRV) == FALSE)}
  assert_that(all(names(fixed_extra) %in% extra_param_names))
  inv_365 = 1.0/365.0
  pts_year = 365.0/time_inc
  n_t_pts = n_years*pts_year
  region = ""

  #Initial SEIRV values
  P0 = array(0, dim = c(n_regions, N_age))
  for(i in 1:n_regions){for(j in 1:N_age){ P0[i,j] = max(1.0, pop_data[i,1, j]) }}

  if(mode_start ==  2){ #Use supplied SEIRV data
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    E_0 = I_0 = array(0, dim = c(n_regions, N_age))
    vacc_initial = array(vacc_data[, 1, ],dim = c(n_regions,N_age))
    V_0 = P0*vacc_initial
    if(mode_start ==  0){ #No initial immunity
      S_0 = P0*(1.0 - vacc_initial)
      R_0 = E_0
    } #if mode_start == 1, S_0 and R_0 set based on FOI/R0 - see packer_setup
  }

  #Population/vaccination
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(n_regions, N_age, n_years))
  for(n_year in 1:n_years){
    for(i in 1:n_regions){
      for(j in 1:N_age){
        dP1_all[i,j,n_year] = max(1.0, pop_data[i, n_year + 1, j])*inv_365
        dP2_all[i,j,n_year] = max(1.0, pop_data[i, n_year, j])*inv_365
        if(j ==  1){
          vacc_rates[i,j,n_year] = vacc_data[i, n_year + 1, j]*inv_365
        } else {
          vacc_rates[i,j,n_year] = max(0.0, vacc_data[i, n_year + 1, j] - vacc_data[i, n_year, j - 1])*inv_365
        }
      }
    }
  }

  #Region and output grouping
  sero_region_groups = sort(unique(sero_template$region))
  sero_regions = regions_breakdown(sero_region_groups)
  case_region_groups = sort(unique(case_template$region))
  case_regions = regions_breakdown(case_region_groups)
  regions_all = regions_breakdown(c(sero_template$region,case_template$region))
  n_regions = length(regions_all)
  i_sero_regions = i_case_regions = rep(0,n_regions)
  i_sero_regions[regions_all %in% sero_regions] = 1
  i_case_regions[regions_all %in% case_regions] = 1

  n_lines_sero = nrow(sero_template)
  unique_lines_sero = unique_lines_case = c(1)
  for(i in 2:n_lines_sero){
    unique = FALSE
    if(sero_template$region[i] %in% sero_template$region[unique_lines_sero]){
      if(sero_template$age_min[i] %in% subset(sero_template[c(1:(i-1)),],
                                              region == sero_template$region[i])$age_min ==  FALSE ||
         sero_template$age_max[i] %in% subset(sero_template[c(1:(i-1)),],
                                              region == sero_template$region[i])$age_max ==  FALSE){unique = TRUE}
    } else {
      unique = TRUE
    }
    if(unique){unique_lines_sero = append(unique_lines_sero,i,after = length(unique_lines_sero))}
  }
  n_sero_pts = length(unique_lines_sero)
  n_lines_case = nrow(case_template)
  for(i in 2:n_lines_case){
    unique = FALSE
    if(case_template$region[i] %in% case_template$region[unique_lines_case] == FALSE){unique = TRUE}
    if(unique){unique_lines_case = append(unique_lines_case,i,after = length(unique_lines_case))}
  }
  n_case_pts = length(unique_lines_case)
  region_index_sero = array(0,dim = c(n_sero_pts,n_regions))
  region_index_case = array(0,dim = c(n_case_pts,n_regions))
  sero_i_age_min = sero_i_age_max = sero_vc_factor = rep(NA,n_sero_pts)
  for(i in 1:n_sero_pts){
    regions_sero_pt = regions_breakdown(sero_template$region[unique_lines_sero[i]])
    for(j in 1:n_regions){
      if(regions_all[j] %in% regions_sero_pt){region_index_sero[i,j] = 1}
    }
    sero_i_age_min[i] = sero_template$age_min[unique_lines_sero[i]]+1
    sero_i_age_max[i] = sero_template$age_max[unique_lines_sero[i]]+1
    sero_vc_factor[i] = sero_template$vc_factor[unique_lines_sero[i]]
  }
  for(i in 1:n_case_pts){
    regions_case_pt = regions_breakdown(case_template$region[unique_lines_case[i]])
    for(j in 1:n_regions){
      if(regions_all[j] %in% regions_case_pt){region_index_case[i,j] = 1}
    }
  }
  output = list(n_regions = n_regions, N_age = N_age, mode_start = mode_start,
                n_years = n_years, n_t_pts = n_t_pts, year0 = year0, time_inc = time_inc,
                t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious,
                E_0 = E_0, I_0 = I_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
                vacc_rate_daily = vacc_rates,
                n_sero_pts = n_sero_pts,n_case_pts = n_case_pts, sero_vc_factor = sero_vc_factor,
                sia_min = sero_i_age_min, sia_max = sero_i_age_max,
                region_index_sero = region_index_sero, region_index_case = region_index_case,
                sero_regions = i_sero_regions,case_regions = i_case_regions, ref_BRA = ref_BRA,
                i_ttf = i_ttf, i_ptf = i_ptf)
  if(mode_start %in% c(0,2)){
    output$S_0 = S_0
    output$R_0 = R_0
  }else{
    output$P0 = P0
  }
  for(name in extra_param_names){
    if(is.null(fixed_extra[[name]]) == FALSE){output[[name]] = fixed_extra[[name]]}
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title fit_data_setup
#'
#' @description Create fitting dataset from serological and case datasets
#'
#' @details TBA
#'
#' @param sero_template Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param case_template Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param year0 First year in population/vaccination data
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param region_index_sero TBA
#' @param region_index_case TBA
#'
#' @export
#'
fit_data_setup <- function(sero_template = data.frame(),case_template = data.frame(), year0 = 1940, time_inc = 5.0,
                           region_index_sero = list(), region_index_case = list()){

  assert_that(is.data.frame(sero_template))
  assert_that(is.data.frame(case_template))
  assert_that(all(c("region","year","age_min","age_max","positives","samples","vc_factor") %in% colnames(sero_template)))
  assert_that(all(c("region","year","cases","deaths") %in% colnames(case_template)))

  year = region = age_min = age_max = 0
  years_data <- sort(unique(c(sero_template$year,case_template$year)))
  i_year_begin = years_data[1] - year0 + 1
  i_year_end = max(years_data) + 1 - year0
  time_pts = (years_data-year0+1)*(365/time_inc)
  sero_region_groups = sort(unique(sero_template$region))
  case_region_groups = sort(unique(case_template$region))
  regions_all = regions_breakdown(c(sero_template$region,case_template$region))
  n_regions = length(regions_all)
  assert_that(n_regions == dim(region_index_sero)[2])
  assert_that(n_regions == dim(region_index_case)[2])
  n_sero_pts = dim(region_index_sero)[1]
  n_case_pts = dim(region_index_case)[1]

  n_t_pts = length(time_pts)
  sero_data_list1 = sero_data_list2 = case_data_list1 = case_data_list2 = list()
  for(i in 1:n_t_pts){
    sero_data_list1[[i]] = sero_data_list2[[i]] = as.numeric(rep(NA,n_sero_pts))
    case_data_list1[[i]] = case_data_list2[[i]] = as.numeric(rep(NA,n_case_pts))
    year_pt = years_data[i]
    sero_subset = subset(sero_template,year == year_pt)
    if(nrow(sero_subset)>0){
      for(j in 1:n_sero_pts){
        regions = regions_all[which(region_index_sero[j,] == 1)]
        region_group = paste(regions,collapse = ",")
        sero_subset2 = subset(sero_subset,region == region_group)
        sero_subset2 = subset(sero_subset2,age_min == sero_template$age_min[j])
        sero_subset2 = subset(sero_subset2,age_max == sero_template$age_max[j])
        sero_data_list1[[i]][j] = sero_subset2$positives[1]
        sero_data_list2[[i]][j] = sero_subset2$samples[1]
      }
    }
    case_subset = subset(case_template,year == year_pt)
    if(nrow(case_subset)>0){
      for(j in 1:n_case_pts){
        regions = regions_all[which(region_index_case[j,] == 1)]
        region_group = paste(regions,collapse = ",")
        k = which(case_subset$region == region_group)
        case_data_list1[[i]][j] = case_subset$cases[k]
        case_data_list2[[i]][j] = case_subset$deaths[k]
      }
    }
  }
  fit_data = data.frame(time = time_pts,obs_sero_positives = I(sero_data_list1), obs_sero_samples = I(sero_data_list2),
                        obs_case_values = I(case_data_list1), obs_death_values = I(case_data_list2))

  return(fit_data)
}
#-------------------------------------------------------------------------------
#' @title packer_setup
#'
#' @description Create packer for fitting
#'
#' @details TBA
#'
#' @param pars_fixed Fixed parameters list created using pars_fixed_setup()
#' @param env_covar_values Environmental covariate values (TBC)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' @param vars_extra_names Names of additional varied parameters
#'
#' @export
#'
packer_setup <- function(pars_fixed = list(), env_covar_values = list(),
                         mode_time = 0, vars_extra_names = c("")){

  #TODO - add assert_that functions?
  #TODO - remove checks/calcs to be moved to FOI/R0 calculation function
  assert_that(all(vars_extra_names %in% extra_param_names))
  if("m_FOI_BRA" %in% vars_extra_names){flag_BRA = 2} else {
    if(is.null(pars_fixed$m_FOI_BRA) == FALSE){ flag_BRA = 1 } else {flag_BRA = 0}
  }
  if(flag_BRA>0){assert_that(is.null(pars_fixed$ref_BRA) == FALSE)}
  assert_that(length(dim(env_covar_values)) == 3)
  n_env_vars = dim(env_covar_values)[1]
  assert_that(dim(env_covar_values)[1] == n_env_vars)
  assert_that(mode_time %in% c(0:5), msg = "mode_time must be an integer between 0 and 5")
  assert_that(all(env_covar_values >=  0.0))
  n_regions = pars_fixed$n_r
  assert_that(dim(env_covar_values)[2] == n_regions)
  time_inc = pars_fixed$time_inc
  pts_year = 365.0/time_inc
  n_years = pars_fixed$n_years
  n_t_pts = n_years*pts_year
  inv_365 = 1.0/365.0
  n_req = switch(mode_time + 1, 1, n_years, 12, pts_year, n_years*12, n_t_pts)
  assert_that(dim(env_covar_values)[3] == n_req)
  for(name in extra_param_names_req){
    if(is.null(pars_fixed[[name]])){assert_that(name %in% vars_extra_names)}
  }
  if(is.null(pars_fixed$i_ttf)==FALSE){
    for(name in extra_param_names_ttf){
      if(is.null(pars_fixed[[name]])){assert_that(name %in% vars_extra_names)}
    }
  }
  if(is.null(pars_fixed$ref_BRA)==FALSE && is.null(pars_fixed$m_FOI_BRA)){
    assert_that("m_FOI_BRA" %in% vars_extra_names)}
  if(is.null(pars_fixed$i_ptf)==FALSE && is.null(pars_fixed$log_a_ptf)){
    assert_that("log_a_ptf" %in% vars_extra_names)}

  date_values = switch(mode_time + 1,
                       rep(1, n_t_pts),
                       sort(rep(c(1:n_years), pts_year)),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12),
                       1 + (floor(time_inc*c(0:(n_t_pts - 1))) %% pts_year),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12) +
                         (12*sort(rep(c(1:n_years) - 1,pts_year))),
                       c(1:n_t_pts))

  packer <- monty_packer(scalar = vars_extra_names,
                         array = list("log_FOI_coeffs" = n_env_vars,
                                      "log_R0_coeffs" = n_env_vars),
                         fixed = pars_fixed,
                         process = function(p){

                           epi_params = epi_param_calc2(pars_fixed, env_covar_values, p$log_FOI_coeffs,
                                                        p$log_R0_coeffs,p)
                           FOI_spillover_t = epi_params$FOI_spillover[,date_values]
                           R0_t = epi_params$R0[,date_values]

                           params_list = list(FOI_spillover = FOI_spillover_t,R0 = R0_t)
                           if(pars_fixed$mode_start == 1){
                             vacc_initial = pars_fixed$V_0
                             N_age = pars_fixed$N_age
                             P0 = pars_fixed$P0
                             S_0 = R_0 = array(0, dim = c(n_regions, N_age))
                             ages = c(1:N_age) - 1
                             for(i in 1:n_regions){
                               R0_year0 = mean(R0_t[i,c(1:pts_year)])
                               FOI_spillover_year0 = mean(FOI_spillover_t[i,c(1:pts_year)])
                               if(R0_year0 <=  1.0){
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
                                   R_0[i,j] = 0
                                   S_0[i,j] = P0[i,j]*(1.0 - vacc_initial[i,j])
                                 }
                               }
                             }
                             params_list$S_0 = S_0
                             params_list$R_0 = R_0
                           }
                           return(params_list)
                         })
  return(packer)
}
#-------------------------------------------------------------------------------
#' @title prior_setup
#'
#' @description Set up prior function for fitting
#'
#' @details TBA
#'
#' @param packer Packer created using packer_setup()
#' @param env_covar_values Environmental covariate values (TBC)
#' @param pars_var Data frame of information on varied parameters, created using pars_var_setup()
#' @param FOI_R0_prior_data Data frame of prior data for FOI/R0
#'
#' @export
#'
prior_setup <- function(packer = NULL, env_covar_values = list(), pars_var = list(), FOI_R0_prior_data = list()){

  #TODO - assert_that functions
  n_env_vars = dim(env_covar_values)[1]
  n_regions = dim(env_covar_values)[2]
  n_req = dim(env_covar_values)[3]
  pts1 = which(grepl("log_FOI_coeffs",pars_var$name))
  pts2 = which(grepl("log_R0_coeffs",pars_var$name))
  pts_add = which(c(1:nrow(pars_var)) %in% c(pts1,pts2) ==  FALSE)
  vars_extra_names = pars_var$name[pts_add]
  assert_that(all(vars_extra_names %in% extra_param_names))
  n_add = length(pts_add)
  env_covars_mean = array(NA,dim = c(n_env_vars,n_regions))
  for(i in 1:n_regions){
    env_covars_mean[,i] = rowMeans(array(env_covar_values[,i,],dim = c(n_env_vars,n_req)))
  }

  prior_function <- function(log_FOI_coeffs = rep(-15,5),log_R0_coeffs = rep(-5,5),
                             p_rep_severe = 0,p_rep_death = 0,p_severe_inf = 0.12,p_death_severe_inf = 0.39,
                             vaccine_efficacy = 1.0,m_FOI_BRA = 1.0,log_a_ptf=1.0,
                             a_T0 = 2.248952, a_Tm = 40.13383, a_c = 0.000271964,
                             mu_T0 = 12.71508, mu_Tm = 38.04809469, mu_c = -0.757869,
                             PDR_T0 = 17.33263, PDR_Tm = 42.19592, PDR_c = 0.000135891,...){

    #Prior applied to coefficients of environmental covariates
    prior_lfc = prior_lrc = rep(0,n_env_vars)
    for(i in 1:n_env_vars){
      j = pts1[i]
      k = pts2[i]
      prior_lfc = log(dtrunc(x = log_FOI_coeffs[i],spec = "norm",a = pars_var$min[j],b = pars_var$max[j],
                             mean = pars_var$mean[j],sd = pars_var$sd[j]))
      prior_lrc = log(dtrunc(x = log_R0_coeffs[i],spec = "norm",a = pars_var$min[k],b = pars_var$max[k],
                             mean = pars_var$mean[k],sd = pars_var$sd[k]))
    }

    #Prior applied to mean FOI and R0
    FOI_mean = R0_mean = rep(NA,n_regions)
    for(i in 1:n_regions){
      FOI_mean[i] = sum(env_covars_mean[,i]*exp(log_FOI_coeffs))
      R0_mean[i] = sum(env_covars_mean[,i]*exp(log_R0_coeffs))
    }
    prior_FOI = log(dtrunc(x = FOI_mean,spec = "norm",a = FOI_R0_prior_data$min[1],b = FOI_R0_prior_data$max[1],
                           mean = FOI_R0_prior_data$mean[1],sd = FOI_R0_prior_data$sd[1]))
    prior_R0 = log(dtrunc(x = R0_mean,spec = "norm",a = FOI_R0_prior_data$min[2],b = FOI_R0_prior_data$max[2],
                          mean = FOI_R0_prior_data$mean[2],sd = FOI_R0_prior_data$sd[2]))

    #Prior applied to additional parameters
    prior_add = rep(0,n_add)
    for(i in 1:n_add){
      j = pts_add[i]
      prior_add[i] = log(dtrunc(x = get(vars_extra_names[i]),spec = "norm",a = pars_var$min[j],
                                b = pars_var$max[j],mean = pars_var$mean[j],sd = pars_var$sd[j]))
    }

    return(sum(c(prior_lfc,prior_lrc,prior_FOI,prior_R0,prior_add)))
  }
  prior <- monty_model_function(prior_function,packer = packer)

  return(prior)
}
#-------------------------------------------------------------------------------
#' @title m_prelim_fit
#'
#' @description Function for preliminary estimation of parameter values
#'
#' @details TBA
#'
#' @param fit_data Fitting data created using fit_data_setup()
#' @param packer Packer created using packer_setup()
#' @param prior Prior function created using prior_setup()
#' @param FOI_R0_prior_data Data frame of prior data for FOI/R0
#' @param pars_var Data frame of information on varied parameters, created using pars_var_setup()
#' @param env_covar_values Environmental covariate values (TBC)
#' @param n_steps Number of times to run estimation cycle
#' @param n_iterations Number of iterations to run per cycle
#' @param n_bounds Number of iterations (ones giving highest posterior likelihood) to use to establish bounds
#'  for next cycle
#' @param deterministic TBA
#' @param n_particles Number of particles
#' @param n_threads Number of threads
#' @param seed Random seed (set to NULL if unused)
#'
#' @export
#'
m_prelim_fit <- function(fit_data = list(), packer = NULL, prior = NULL, FOI_R0_prior_data = list(),
                         pars_var = list(), env_covar_values = list(), n_steps = 1, n_iterations = 10,
                         n_bounds = 10, deterministic = FALSE, n_particles = 1, n_threads = 1, seed = NULL){

  #TODO - Add assert_that checks
  #check fit_data
  #check packer
  #check priors
  #check pars_var
  #check env_covar_values
  assert_that(is.logical(deterministic))

  if(deterministic){
    assert_that(n_particles==1 && n_threads==1)
    filter <- dust_unfilter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                                   n_particles = n_particles, n_threads = n_threads)
  } else {
    filter <- dust_filter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                                 n_particles = n_particles, n_threads = n_threads, seed = seed)
  }
  likelihood <- dust_likelihood_monty(obj = filter, packer = packer)
  if(is.null(prior)){posterior <- likelihood}else{posterior<-likelihood+prior}

  n_env_vars = dim(env_covar_values)[1]
  n_regions = dim(env_covar_values)[2]
  n_req = dim(env_covar_values)[3]
  n_params = nrow(pars_var)
  pts1 = which(grepl("log_FOI_coeffs",pars_var$name))
  pts2 = which(grepl("log_R0_coeffs",pars_var$name))

  explore = list()
  for(step in 1:n_steps){
    if(step==1){
      pars_min = pars_var$min
      pars_max = pars_var$max
    } else{
      pars_min = rowMins(t(explore[[step-1]][1:n_bounds, c(1:n_params)]))
      pars_max = rowMaxs(t(explore[[step-1]][1:n_bounds, c(1:n_params)]))
    }
    cat("\n")
    set.seed(seed)
    param_sets = lhs(n_iterations,rect = matrix(c(pars_min,pars_max),ncol = 2))
    explore[[step]] = data.frame(array(NA,dim = c(n_iterations,n_params+1)))
    explore[[step]][,c(1:n_params)] = param_sets
    #TODO - put correct parameter names as headings
    colnames(explore[[step]]) = c(paste0("param",c(1:n_params)),"density")
    for(iter in 1:n_iterations){
      cat("\n",step,"-",iter,":\n",sep="")
      cat(signif(param_sets[iter,],3))
      FOI_mean = R0_mean = rep(NA,n_regions)
      for(n_region in 1:n_regions){
        env_covars_mean = rowMeans(array(env_covar_values[,n_region,],dim = c(n_env_vars,n_req)))
        FOI_mean[n_region] = sum(env_covars_mean*exp(param_sets[iter,pts1]))
        R0_mean[n_region] = sum(env_covars_mean*exp(param_sets[iter,pts2]))
      }
      test1 = any(FOI_mean>FOI_R0_prior_data$max[1])
      test2 = any(R0_mean>FOI_R0_prior_data$max[2])
      if(any(test1,test2)){
        explore[[step]]$density[iter] <- -Inf
        cat("\nRejected (outwith FOI/R0 range)")
      } else {
        explore[[step]]$density[iter] <- monty_model_density(posterior, parameters = param_sets[iter,])
        cat("\nDensity: ",format(signif(explore[[step]]$density[iter], 3), scientific=TRUE))
      }
    }
    cat("\n")
    explore[[step]] = explore[[step]][order(-explore[[step]]$density),]
  }

  return(explore)
}
#-------------------------------------------------------------------------------
#' @title m_sample
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param fit_data Fitting data created using fit_data_setup()
#' @param packer Packer created using packer_setup()
#' @param prior Prior function created using prior_setup()
#' @param FOI_R0_prior_data Data frame of prior data for FOI/R0
#' @param pars_var Data frame of information on varied parameters, created using pars_var_setup()
#' @param initial TBA
#' @param vcv Variance-covariance matrix
#' @param rerun_every TBA
#' @param n_chains Number of chains
#' @param n_iterations Number of iterations for which to run each chain
#' @param deterministic TBA
#' @param n_particles Number of particles
#' @param n_threads Number of threads
#' @param parallel TBA
#' @param output_file Name of file location to save results (set to NULL if unused)
#' @param seed Random seed (set to NULL if unused)
#'
#' @export
#'
m_sample <- function(fit_data = list(), packer = NULL, prior = NULL, FOI_R0_prior_data = list(),
                     pars_var = list(), initial = list(), vcv = list(), rerun_every = 50, n_chains = 1,
                     n_iterations = 10, deterministic = FALSE, n_particles = 1, n_threads = 1,
                     parallel = FALSE, output_file = NULL, seed = NULL){

  #TODO - add assert_that functions
  n_params = nrow(pars_var)
  assert_that(all(dim(vcv)==c(n_params,n_params)))

  if(deterministic){
    assert_that(n_particles==1)
    assert_that(n_threads==1)
    assert_that(rerun_every==Inf,msg="rerun_every must be Inf when in deterministic mode")
    filter <- dust_unfilter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                                 n_particles = n_particles, n_threads = n_threads)
  } else {
    filter <- dust_filter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                                 n_particles = n_particles, n_threads = n_threads, seed = seed)
  }

  likelihood <- dust_likelihood_monty(obj = filter, packer = packer)
  if(is.null(prior)){posterior <- likelihood}else{posterior <- likelihood + prior}
  posterior$domain[,1] = pars_var$min
  posterior$domain[,2] = pars_var$max
  assert_that(all(initial>= posterior$domain[,1]))
  assert_that(all(initial<= posterior$domain[,2]))
  sampler <- monty_sampler_random_walk(vcv = vcv, boundaries = "reject",
                                       rerun_every = rerun_every, rerun_random = FALSE)

  if(parallel){
    prior0 = prior$density(initial[,1])#Seems to be necessary for some reason
    runner = monty_runner_callr(n_workers = n_chains)
  }else{
    runner = monty_runner_serial()}

  set.seed(seed)
  samples <- monty_sample(model = posterior,sampler = sampler, n_steps = n_iterations,
                          initial = initial, n_chains = n_chains, runner = runner)

  if(is.null(output_file)==FALSE){
    saveRDS(list(samples = samples,
                 params = list(pars_var = pars_var, n_chains = n_chains,
                               n_iterations = n_iterations,
                               n_particles = n_particles, n_threads = n_threads)),
            file = output_file)
  }

  return(samples)
}
