#TODO - adapt to allow for case and sero data from same region
#' @title Generate_Dataset
#'
#' @description Generate dataset [TBA]
#'
#' @details This function is used to [TBA]
#'
#' [TBA - Explanation of breakdown of regions to model and how to set lengths of FOI_values and R0_values]
#'
#' @param FOI_values Array of values of force of infection due to spillover from sylvatic reservoir by region + time point
#' @param R0_values Array of values of basic reproduction number for human-human transmission by region and time point
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param template TBA \cr
#' sero: Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#' and number of samples \cr
#' case: Annual reported case/death data template - data frame with region and year \cr
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
#' @param seed Optional random seed value; set to NULL to omit.
#' @param region_grouping TBA
#' @param mode_grouping TBA
#' '
#' @export
#'
Generate_Dataset <- function(FOI_values = c(),R0_values = c(),input_data = list(),
                             template = list(sero=NULL,case=NULL,xref_sero=NULL,xref_case=NULL),
                             vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = NULL, mode_time = 0,
                             n_reps = 1,deterministic = FALSE, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                             p_rep_severe = 1.0,p_rep_death = 1.0,mode_parallel = FALSE,cluster = NULL,output_frame = FALSE,
                             seed = NULL, region_grouping=NULL, mode_grouping=1){

  assert_that(input_data_check(input_data),msg = paste("Input data must be in standard format",
                                                       " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  #TBA - Change assert_that functions for template
  assert_that(any(is.null(template$sero) == FALSE,is.null(template$case) == FALSE),
              msg = "Need serological and/or case data template(s)")
  if(is.null(template$sero) == FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(template$sero)))
  }
  if(is.null(template$case) == FALSE){ #TODO - add option for monthly data
    assert_that(all(c("region","year","cases","deaths") %in% names(template$case)))
    #if("month" %in% names(template$case)){flag_monthly_cases=TRUE}else{flag_monthly_cases=FALSE}
    assert_that(between(p_severe_inf,0.0,1.0),msg = "Severe infection rate must be between 0-1")
    assert_that(between(p_death_severe_inf,0.0,1.0),msg = "Fatality rate of severe infections must be between 0-1")
    assert_that(between(p_rep_severe,0.0,1.0),msg = "Severe infection reporting probability must be between 0-1")
    assert_that(between(p_rep_death,0.0,1.0),msg = "Fatal infection reporting probability must be between 0-1")
  }
  assert_that(between(vaccine_efficacy,0.0,1.0),msg = "Vaccine efficacy must be between 0-1")
  assert_that(is.logical(mode_parallel))
  if(mode_parallel){assert_that(is.null(cluster) == FALSE)}
  assert_that(length(dim(FOI_values)) == 2,msg = "FOI_values must be 2-D array")
  assert_that(length(dim(R0_values)) == 2,msg = "R0_values must be 2-D array")

  #Check that regions in template match up with input data
  regions = regions_breakdown(c(template$sero$region,template$case$region))
  assert_that(all(regions==input_data$region_labels),
              msg="Regions in input data must match regions in template; use input_data_truncate() to adjust input data")
  n_regions = length(regions)
  assert_that(dim(FOI_values)[1] == n_regions,msg = "1st dimension of FOI_values must match number of regions to be modelled")
  assert_that(dim(R0_values)[1] == n_regions,msg = "1st dimension of R0_values must match number of regions to be modelled")
  if(mode_start == 2){assert_that(length(start_SEIRV) == n_regions,
                                  msg = "Number of start_SEIRV datasets must match number of regions")}

  #Group regions based on template
  if(is.null(region_grouping)){
    region_grouping = get_region_grouping(regions,template,mode_grouping)
  }
  sero_line_list=region_grouping$sero_line_list
  case_line_list=region_grouping$case_line_list
  n_groups=length(region_grouping$region_groups)

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(template$sero)){model_sero_data = NULL} else {
    blank = rep(0,nrow(template$sero))
    model_sero_data = data.frame(samples = blank,positives = blank,sero = blank)
  }
  if(is.null(template$case)){
    model_case_values = model_death_values = NA
  } else {
    model_case_values = model_death_values = rep(0,nrow(template$case))
  }

  n_years=length(input_data$years_labels)
  N_age=length(input_data$age_labels)
  if(mode_parallel){
    FOI_subsets = R0_subsets = vacc_data_subsets = pop_data_subsets = years_data_sets = start_SEIRV_sets = list()
    for(n_group in 1:n_groups){
      i_regions=region_grouping$region_groups[[n_group]]
      n_regions2=length(i_regions)
      FOI_subsets[[n_group]] = array(FOI_values[i_regions,],dim=c(n_regions2,dim(FOI_values)[2]))
      R0_subsets[[n_group]] = array(R0_values[i_regions,],dim=c(n_regions2,dim(R0_values)[2]))
      vacc_data_subsets[[n_group]] = array(input_data$vacc_data[i_regions,,],dim=c(n_regions2,n_years,N_age))
      pop_data_subsets[[n_group]] = array(input_data$pop_data[i_regions,,],dim=c(n_regions2,n_years,N_age))
      start_SEIRV_sets[[n_group]] = list() #TBA
    }

    model_output_all = clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_subsets, R0 = R0_subsets,
                                  vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                  years_data = region_grouping$years_data, start_SEIRV = start_SEIRV_sets,
                                  mode_out = region_grouping$mode_out,
                                  MoreArgs = list(year0 = input_data$years_labels[1],vaccine_efficacy = vaccine_efficacy,
                                                  time_inc = time_inc,mode_start = mode_start,mode_time = mode_time,
                                                  n_particles = n_reps, n_threads = 1 ,deterministic = deterministic,
                                                  seed = seed))
  }

  for(n_group in 1:n_groups){
    i_regions=region_grouping$region_groups[[n_group]]
    n_regions2=length(i_regions)
    if(n_regions2>0){
      if(mode_parallel){
        model_output = model_output_all[[n_group]]
      } else {
        model_output = Model_Run(FOI_spillover = array(FOI_values[i_regions,],dim=c(n_regions2,dim(FOI_values)[2])),
                                 R0 = array(R0_values[i_regions,],dim=c(n_regions2,dim(R0_values)[2])),
                                 vacc_data = array(input_data$vacc_data[i_regions,,],dim=c(n_regions2,n_years,N_age)),
                                 pop_data = array(input_data$pop_data[i_regions,,],dim=c(n_regions2,n_years,N_age)),
                                 years_data = region_grouping$years_data[[n_group]],
                                 year0 = input_data$years_labels[1], vaccine_efficacy = vaccine_efficacy,
                                 time_inc = time_inc, mode_out = region_grouping$mode_out[[n_group]],mode_start = mode_start,
                                 start_SEIRV = NULL, #TBC
                                 mode_time = mode_time,n_particles = n_reps,
                                 n_threads = n_reps, deterministic = deterministic, seed = seed)
      }

      for(n_region in i_regions){
        n_region2=match(regions[n_region],regions[region_grouping$region_groups[[n_group]]])

        #Compile case data if needed
        if(is.na(case_line_list[[n_region]][1]) == FALSE){
          case_line_list_region = case_line_list[[n_region]]
          years_case = template$case$year[case_line_list_region]
          n_lines = length(case_line_list_region)

          for(n_rep in 1:n_reps){
            rep_cases = rep_deaths = rep(0,n_lines)
            for(n_line in 1:n_lines){
              #pts = c(1:t_pts)[model_output$year == years_case[n_line]]
              pts = which(model_output$year==years_case[n_line])
              infs = model_output$C_annual[n_region2,n_rep,pts]
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
            sero_results = sero_calculate2_alt(template$sero[sero_line_list_region,],
                                               model_output, n_region2, n_rep)
            model_sero_data$samples[sero_line_list_region] = model_sero_data$samples[sero_line_list_region]+sero_results$samples
            model_sero_data$positives[sero_line_list_region] = model_sero_data$positives[sero_line_list_region] + sero_results$positives
          }
        }
      }
    }
  }

  if(is.null(template$sero) == FALSE){model_sero_data$sero = model_sero_data$positives/model_sero_data$samples}
  if(is.null(template$case) == FALSE && n_reps>1){
    model_case_values = round(model_case_values/n_reps)
    model_death_values = round(model_death_values/n_reps)
  }

  output = list()
  if(output_frame) { #Output complete frames of data
    if(is.null(template$sero) == FALSE){
      output$model_sero_data = data.frame(region = template$sero$region,
                                          year = template$sero$year,
                                          age_min = template$sero$age_min,
                                          age_max = template$sero$age_max,
                                          samples = template$sero$samples,
                                          positives = template$sero$samples*model_sero_data$sero,
                                          vc_factor = template$sero$vc_factor)

    }
   if(is.null(template$case) == FALSE){
     output$model_case_data = data.frame(region = template$case$region,
                                         year = template$case$year,
                                         cases = model_case_values,
                                         deaths = model_death_values)

   }
  } else { #Minimal output
    output$model_sero_values = model_sero_data$sero
    output$model_case_values = model_case_values
    output$model_death_values = model_death_values
  }

  return(output)
}
