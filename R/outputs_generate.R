# R file for functions used to generate sets of annual serological data, annual case/death data
# and annual burden (VIMC format) data
#-------------------------------------------------------------------------------
#' @title Generate_Dataset
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
#' '
#' @export
#'
Generate_Dataset <- function(FOI_values = c(),R0_values = c(),input_data = list(),sero_template = NULL,case_template = NULL,
                             vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = NULL, mode_time = 0,
                             n_reps = 1,deterministic = FALSE, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                             p_rep_severe = 1.0,p_rep_death = 1.0,mode_parallel = FALSE,cluster = NULL,output_frame = FALSE,
                             seed = NULL){

  assert_that(input_data_check(input_data),msg = paste("Input data must be in standard format",
                                                     " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template) == FALSE,is.null(case_template) == FALSE),msg = "Need at least one template")
  if(is.null(sero_template) == FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template) == FALSE){
    assert_that(all(c("region","year") %in% names(case_template)))
    assert_that(p_severe_inf >= 0.0 && p_severe_inf <= 1.0,msg = "Severe infection rate must be between 0-1")
    assert_that(p_death_severe_inf >= 0.0 && p_death_severe_inf <= 1.0,
                msg = "Fatality rate of severe infections must be between 0-1")
    assert_that(p_rep_severe >= 0.0 && p_rep_severe <= 1.0,msg = "Severe infection reporting probability must be between 0-1")
    assert_that(p_rep_death >= 0.0 && p_rep_death <= 1.0,msg = "Fatal infection reporting probability must be between 0-1")
  }
  assert_that(is.logical(mode_parallel))
  if(mode_parallel){assert_that(is.null(cluster) == FALSE)}

  #Prune input data based on regions
  regions = regions_breakdown(c(sero_template$region,case_template$region))
  input_data = input_data_truncate(input_data,regions)
  n_regions = length(input_data$region_labels)

  #Cross-reference templates with input regions
  if(is.null(sero_template) == FALSE){
    xref_sero = template_region_xref(sero_template,input_data$region_labels)
    sero_line_list = xref_sero$line_list
  } else {
    xref_sero = data.frame(year_data_begin = rep(Inf,n_regions),year_end = rep(-Inf,n_regions))
    sero_line_list = rep(NA,n_regions)
  }
  if(is.null(case_template) == FALSE){
    xref_case = template_region_xref(case_template,input_data$region_labels)
    case_line_list = xref_case$line_list
  } else {
    xref_case = data.frame(year_data_begin = rep(Inf,n_regions),year_end = rep(-Inf,n_regions))
    case_line_list = rep(NA,n_regions)
  }
  year_data_begin = year_end = rep(NA,length(input_data$region_labels))
  for(i in 1:length(year_data_begin)){
    year_data_begin[i] = min(xref_sero$year_data_begin[i],xref_case$year_data_begin[i])
    year_end[i] = max(xref_sero$year_end[i],xref_case$year_end[i])
  }

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
  output_types = rep(NA,n_regions)
  for(n_region in 1:n_regions){
    if(is.na(case_line_list[[n_region]][1]) == FALSE){
      if(is.na(sero_line_list[[n_region]][1]) == FALSE){ output_types[n_region] = "case_sero"}else{output_types[n_region] = "case"}
    } else {
      output_types[n_region] = "sero"
    }
  }

  if(mode_parallel){
    FOI_subsets = R0_subsets = vacc_data_subsets = pop_data_subsets = years_data_sets = list() #TODO - change input data?
    for(n_region in 1:n_regions){
      FOI_subsets[[n_region]] = FOI_values[n_region,]
      R0_subsets[[n_region]] = R0_values[n_region,]
      vacc_data_subsets[[n_region]] = input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]] = input_data$pop_data[n_region,,]
      years_data_sets[[n_region]] = c(year_data_begin[n_region]:year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV = rep(NA,n_regions)}
    #TODO - Test
    model_output_all = clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_subsets, R0 = R0_subsets,
                                  vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                  years_data = years_data_sets, start_SEIRV = start_SEIRV, output_type = output_types,
                                  MoreArgs = list(year0 = input_data$years_labels[1],vaccine_efficacy = vaccine_efficacy,
                                                  time_inc = time_inc,mode_start = mode_start,mode_time = mode_time,
                                                  n_particles = n_reps, n_threads = 1 ,deterministic = deterministic))
  }

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel == FALSE){
      if(is.null(seed) == FALSE && deterministic == FALSE){set.seed(seed)}
      model_output = Model_Run(FOI_spillover = FOI_values[n_region,],R0 = R0_values[n_region,],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(year_data_begin[n_region]:year_end[n_region]), year0 = input_data$years_labels[1],
                               vaccine_efficacy = vaccine_efficacy, time_inc = time_inc, output_type = output_types[n_region],
                               mode_start = mode_start, start_SEIRV = start_SEIRV[[n_region]],mode_time = mode_time,
                               n_particles = n_reps,n_threads = n_reps, deterministic = deterministic)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts = length(model_output$year)

    #Compile case data if needed
    if(is.na(case_line_list[[n_region]][1]) == FALSE){
      case_line_list_region = case_line_list[[n_region]]
      years_case = case_template$year[case_line_list_region]
      n_lines = length(case_line_list_region)

      for(n_rep in 1:n_reps){
        rep_cases = rep_deaths = rep(0,n_lines)
        for(n_line in 1:n_lines){
          pts = c(1:t_pts)[model_output$year == years_case[n_line]]
          infs = sum(model_output$C[n_rep,pts])
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
        sero_results = sero_calculate2(sero_template[sero_line_list_region,],model_output,n_rep)
        model_sero_data$samples[sero_line_list_region] = model_sero_data$samples[sero_line_list_region]+sero_results$samples
        model_sero_data$positives[sero_line_list_region] = model_sero_data$positives[sero_line_list_region] +
          sero_results$positives
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
#-------------------------------------------------------------------------------
#' @title Generate_VIMC_Burden_Dataset
#'
#' @description Generate annual burden data in format used by VIMC
#'
#' @details This function is used to generate annual burden data in the format used by the Vaccine Impact
#'   Modelling Consortium (VIMC) [TBA]
#'
#' @param FOI_values Array of values of force of infection due to spillover from sylvatic reservoir by region + time point
#' @param R0_values Array of values of basic reproduction number for human-human transmission by region and time point
#' @param input_data List of population and vaccination data for multiple regions
#' @param template Burden data in VIMC format, with regions, years, minimum and maximum age, and life expectancy
#'   for each line
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
#' @param YLD_per_case TBA
#' @param mode_parallel TRUE/FALSE - set model to run in parallel using cluster if TRUE
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param seed Optional random seed value to set before running each region for stochastic normalization; set to NULL
#'   to omit; will not work if mode_parallel is not set to FALSE.
#' '
#' @export
#'
Generate_VIMC_Burden_Dataset <- function(FOI_values = c(),R0_values = c(),input_data = list(),template = NULL,
                                         vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = NULL, mode_time = 0,
                                         n_reps = 1,deterministic = FALSE, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                         p_rep_severe = 1.0,p_rep_death = 1.0,YLD_per_case = 0.006486, mode_parallel = FALSE,
                                         cluster = NULL,seed = NULL){

  assert_that(input_data_check(input_data),msg = paste("Input data must be in standard format",
                                                     " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(all(c("region","year","age_min","age_max","life_exp") %in% names(template)))
  assert_that(vaccine_efficacy >= 0.0 && vaccine_efficacy <= 1.0,msg = "Vaccine efficacy must be between 0-1")
  assert_that(p_severe_inf >= 0.0 && p_severe_inf <= 1.0,msg = "Severe infection rate must be between 0-1")
  assert_that(p_death_severe_inf >= 0.0 && p_death_severe_inf <= 1.0,msg = "Severe infection fatality rate must be between 0-1")
  assert_that(is.logical(mode_parallel))
  if(mode_parallel){assert_that(is.null(cluster) == FALSE)}

  #Prune input data based on regions
  regions = regions_breakdown(template$region)
  input_data = input_data_truncate(input_data,regions)
  n_regions = length(input_data$region_labels)
  assert_that(length(dim(FOI_values)) == 2,msg = "FOI_values must be 2-D array")
  assert_that(length(dim(R0_values)) == 2,msg = "R0_values must be 2-D array")
  assert_that(dim(FOI_values)[1] == n_regions,msg = "1st dimension of FOI_values must match number of regions to be modelled")
  assert_that(dim(R0_values)[1] == n_regions,msg = "1st dimension of R0_values must match number of regions to be modelled")
  if(mode_start == 2){assert_that(length(start_SEIRV) == n_regions,
                                  msg = "Number of start_SEIRV datasets must match number of regions")}

  #Cross-reference templates with input regions
  xref = template_region_xref(template,input_data$region_labels)
  n_lines_total = nrow(template)

  #Set up data structures to take modelled data
  model_case_values = model_death_values = cohort_size = rep(0,n_lines_total)

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel){
    FOI_subsets = R0_subsets = vacc_data_subsets = pop_data_subsets = years_data_sets = list() #TODO - change input data?
    for(n_region in 1:n_regions){
      FOI_subsets[[n_region]] = FOI_values[n_region,]
      R0_subsets[[n_region]] = R0_values[n_region,]
      vacc_data_subsets[[n_region]] = input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]] = input_data$pop_data[n_region,,]
      years_data_sets[[n_region]] = c(xref$year_data_begin[n_region]:xref$year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV = rep(NA,n_regions)}
    model_output_all = clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_subsets, R0 = R0_subsets,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV,
                                MoreArgs = list(output_type = "case_alt", year0 = input_data$years_labels[1],
                                              mode_start = mode_start, vaccine_efficacy = vaccine_efficacy,
                                              time_inc = time_inc, n_particles = n_reps, n_threads = 1 ,
                                              deterministic = deterministic, mode_time = mode_time))
  }

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel == FALSE){
      if(is.null(seed) == FALSE && deterministic == FALSE){set.seed(seed)}
      model_output = Model_Run(FOI_spillover = FOI_values[n_region,],R0 = R0_values[n_region,],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(xref$year_data_begin[n_region]:xref$year_end[n_region]),
                               start_SEIRV = start_SEIRV[[n_region]],output_type = "case_alt",
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, time_inc = time_inc, n_particles = n_reps,
                               n_threads = n_reps, deterministic = deterministic, mode_time = mode_time)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts_all = c(1:length(model_output$year))

    #Compile case data
    case_line_list_region = xref$line_list[[n_region]]
    years_case = template$year[case_line_list_region]
    n_lines = length(case_line_list_region)
    age_pts = list()
    for(n_line in 1:n_lines){
      line = case_line_list_region[n_line]
      n_age_min = findInterval(template$age_min[line],input_data$age_labels)
      n_age_max = findInterval(template$age_max[line]-1,input_data$age_labels)
      age_pts[[n_line]] = c(n_age_min:n_age_max)
      cohort_size[line] = sum(input_data$pop_data[n_region,
                                                input_data$years_labels == years_case[n_line],age_pts[[n_line]]])
    }

    for(n_rep in 1:n_reps){
      cases = deaths = rep(0,n_lines)
      for(n_line in 1:n_lines){
        year = years_case[n_line]
        t_pts = t_pts_all[model_output$year == years_case[n_line]]
        infs = sum(model_output$C[age_pts[[n_line]],n_rep,t_pts])
        if(deterministic){
          cases[n_line] = infs*p_severe_inf
          deaths[n_line] = cases[n_line]*p_death_severe_inf
        } else {
          cases[n_line] = rbinom(1,floor(infs),p_severe_inf)
          deaths[n_line] = rbinom(1,cases[n_line],p_death_severe_inf)
        }
      }
      model_case_values[case_line_list_region] = model_case_values[case_line_list_region]+cases
      model_death_values[case_line_list_region] = model_death_values[case_line_list_region]+deaths
    }
  }

  model_case_values = model_case_values/n_reps
  model_death_values = model_death_values/n_reps
  model_YLL_values = model_death_values*template$life_exp
  model_dalys_values = (model_case_values*YLD_per_case)+model_YLL_values

  return(data.frame(disease = rep("YF",n_lines_total),year = template$year,age = template$age_min,
                    region_label1 = substr(template$region,1,3),region_label2 = template$region,
                    cohort_size = cohort_size,cases = model_case_values,dalys = model_dalys_values,
                    deaths = model_death_values,YLL = model_YLL_values))
}
