# R file for functions used to generate sets of annual serological data, annual case/death data
# and annual burden (VIMC format) data
#-------------------------------------------------------------------------------
#' @title Generate_Dataset
#'
#' @description Generate annual serological and/or case/death data
#'
#' @details This function is used to generate annual serological and/or case/death data based on templates;
#' it is normally used by the single_like_calc() function. The separate Generate_Sero_dataset and
#' Generate_Case_Dataset functions can be used when only one or the other type is required; this function exists
#' to cover instances where serological and case data may need to be generated for the same region.
#'
#' [TBA - Explanation of breakdown of regions to model and how to set lengths of FOI_values and R0_values]
#'
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param FOI_values Vector of values of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Vector of values of the basic reproduction number for human-human transmission
#' @param sero_template Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#'   and number of samples
#' @param case_template Annual reported case/death data template - data frame with region and year
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' @param output_frame Flag indicating whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#' '
#' @export
#'
Generate_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),sero_template = NULL,case_template = NULL,
                             vaccine_efficacy = 1.0, p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = 1.0,
                             p_rep_death = 1.0,mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                             mode_parallel = "none",cluster = NULL,output_frame=FALSE){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template)==FALSE,is.null(case_template)==FALSE),msg="Need at least one template")
  if(is.null(sero_template)==FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template)==FALSE){
    assert_that(all(c("region","year") %in% names(case_template)))
    assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
    assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,
                msg="Fatality rate of severe infections must be between 0-1")
    assert_that(p_rep_severe>=0.0 && p_rep_severe<=1.0,msg="Severe infection reporting probability must be between 0-1")
    assert_that(p_rep_death>=0.0 && p_rep_death<=1.0,msg="Fatal infection reporting probability must be between 0-1")
  }
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  #Prune input data based on regions
  regions=regions_breakdown(c(sero_template$region,case_template$region))
  n_regions=length(input_data$region_labels)
  input_data=input_data_truncate(input_data,regions)

  #Cross-reference templates with input regions
  if(is.null(sero_template)==FALSE){
    xref_sero=template_region_xref(sero_template,input_data$region_labels)
    sero_line_list=xref_sero$line_list
  } else {
    xref_sero=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    sero_line_list=rep(NA,n_regions)
  }
  if(is.null(case_template)==FALSE){
    xref_case=template_region_xref(case_template,input_data$region_labels)
    case_line_list=xref_case$line_list
  } else {
    xref_case=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    case_line_list=rep(NA,n_regions)
  }
  year_data_begin=year_end=rep(NA,length(input_data$region_labels))
  for(i in 1:length(year_data_begin)){
    year_data_begin[i]=min(xref_sero$year_data_begin[i],xref_case$year_data_begin[i])
    year_end[i]=max(xref_sero$year_end[i],xref_case$year_end[i])
  }

  inv_reps=1/n_reps
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions to be modelled")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions to be modelled")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(sero_template)){model_sero_data=NULL} else {
    blank=rep(0,nrow(sero_template))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  }
  if(is.null(case_template)){model_case_values=model_death_values=NA} else {
    model_case_values=model_death_values=rep(0,nrow(case_template))
  }

  #Set up vector of output types to get from model if needed
  if(mode_parallel %in% c("none","clusterMap")){
    output_types=rep(NA,n_regions)
    for(n_region in 1:n_regions){
      if(is.na(case_line_list[[n_region]][1])==FALSE){
        if(is.na(sero_line_list[[n_region]][1])==FALSE){output_types[n_region]="case+sero"} else{output_types[n_region]="case"}
      } else {output_types[n_region]="sero"}
    }
  }

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(year_data_begin):max(year_end))
    if(is.null(sero_template)==FALSE){if(is.null(case_template)==FALSE){output_type="case+sero"} else {output_type="sero"}
      } else {output_type="case"}
    model_output_all=Model_Run_Multi_Input(FOI_spillover = FOI_values,R0 = R0_values,
                                           vacc_data = input_data$vacc_data, pop_data = input_data$pop_data,
                                           years_data = years_data_all, start_SEIRV=start_SEIRV,output_type = output_type,
                                           year0 = input_data$years_labels[1],mode_start = mode_start,
                                           vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                           n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(year_data_begin[n_region]:year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV, output_type = output_types,
                                MoreArgs=list(year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(year_data_begin[n_region]:year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = output_types[n_region],
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
      #cat("\n\t\tFinished modelling region ",n_region)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile case data if needed
    if(is.na(case_line_list[[n_region]][1])==FALSE){
      case_line_list_region=case_line_list[[n_region]]
      years_case=case_template$year[case_line_list_region]
      n_lines=length(case_line_list_region)

      for(n_rep in 1:n_reps){
        rep_cases=rep_deaths=rep(0,n_lines)
        for(n_line in 1:n_lines){
          pts=c(1:t_pts)[model_output$year==years_case[n_line]]
          infs=sum(model_output$C[n_rep,pts])
          if(deterministic){
            severe_infs=floor(infs)*p_severe_inf
            deaths=severe_infs*p_death_severe_inf
            rep_deaths[n_line]=round(deaths*p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)

          } else {
            severe_infs=rbinom(1,floor(infs),p_severe_inf)
            deaths=rbinom(1,severe_infs,p_death_severe_inf)
            rep_deaths[n_line]=rbinom(1,deaths,p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
          }
        }

        model_case_values[case_line_list_region]=model_case_values[case_line_list_region]+rep_cases
        model_death_values[case_line_list_region]=model_death_values[case_line_list_region]+rep_deaths
      }
    }

    #Compile seroprevalence data if necessary
    if(is.na(sero_line_list[[n_region]][1])==FALSE){
      sero_line_list_region=sero_line_list[[n_region]]
      for(n_rep in 1:n_reps){
        sero_results=sero_calculate2(sero_template[sero_line_list_region,],model_output,n_rep)
        model_sero_data$samples[sero_line_list_region]=model_sero_data$samples[sero_line_list_region]+sero_results$samples
        model_sero_data$positives[sero_line_list_region]=model_sero_data$positives[sero_line_list_region]+sero_results$positives
      }
    }
  }

  if(is.null(sero_template)==FALSE){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(is.null(case_template)==FALSE){
    model_case_values=model_case_values*inv_reps
    model_death_values=model_death_values*inv_reps
  }

  if(output_frame) { #Output complete frames of data
    return(list(model_sero_data=data.frame(region=sero_template$region,year=sero_template$year,
                                           age_min=sero_template$age_min,age_max=sero_template$age_max,
                                           samples=sero_template$samples,positives=sero_template$samples*model_sero_data$sero),
                model_case_data=data.frame(region=case_template$region,year=case_template$year,
                                           cases=model_case_values,deaths=model_death_values)))
  } else { #Minimal output for MCMC
    return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
                model_death_values=model_death_values))
  }
}
#-------------------------------------------------------------------------------
#' @title Generate_Sero_Dataset
#'
#' @description Generate serological dataset
#'
#' @details This function is used to generate serological, data for multiple regions/parameter sets based
#'   on a template.
#'
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param template Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' @param output_frame Flag indicating whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#' '
#' @export
#'
Generate_Sero_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),template = NULL,
                                  vaccine_efficacy = 1.0, mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1,
                                  deterministic = FALSE, mode_parallel = "none",cluster = NULL, output_frame=FALSE){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  assert_that(all(c("age_min","age_max","samples","positives","vc_factor","region") %in% names(template)))
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0-1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  #Prune input data based on regions
  regions=regions_breakdown(template$region)
  input_data=input_data_truncate(input_data,regions)

  #Cross-reference templates with input regions
  xref=template_region_xref(template,input_data$region_labels)

  inv_reps=1/n_reps
  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data
  blank=rep(0,nrow(template))
  model_sero_data=list(samples=blank,positives=blank,sero=blank)

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(xref$year_data_begin):max(xref$year_end))
    model_output_all=Model_Run_Multi_Input(FOI_spillover = FOI_values,R0 = R0_values,
                                           vacc_data = input_data$vacc_data, pop_data = input_data$pop_data,
                                           years_data = years_data_all, start_SEIRV=start_SEIRV,output_type = "sero",
                                           year0 = input_data$years_labels[1],mode_start = mode_start,
                                           vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                           n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(xref$year_data_begin[n_region]:xref$year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV,
                                MoreArgs=list(output_type = "sero", year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  # Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(xref$year_data_begin[n_region]:xref$year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = "sero",
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile seroprevalence data
    sero_line_list_region=xref$line_list[[n_region]]
    for(n_rep in 1:n_reps){
      sero_results=sero_calculate2(template[sero_line_list_region,],model_output,n_rep)
      model_sero_data$samples[sero_line_list_region]=model_sero_data$samples[sero_line_list_region]+sero_results$samples
      model_sero_data$positives[sero_line_list_region]=model_sero_data$positives[sero_line_list_region]+sero_results$positives
    }
  }
  model_sero_data$sero=model_sero_data$positives/model_sero_data$samples

  if(output_frame){ #Output complete frame of data
    return(data.frame(region=template$region,year=template$year,
                      age_min=template$age_min,age_max=template$age_max,
                      samples=template$samples,positives=template$samples*model_sero_data$sero))
  } else
  { #Minimal output for MCMC
    return(model_sero_data)
  }
}
#-------------------------------------------------------------------------------
#' @title Generate_Case_Dataset
#'
#' @description Generate annual case/death data
#'
#' @details This function is used to generate annual case/death data based on a template
#'
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param template Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' @param output_frame Flag indicating whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#'
#' @export
#'
Generate_Case_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),template = NULL,vaccine_efficacy = 1.0,
                                  p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = 1.0,p_rep_death = 1.0,
                                  mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                                  mode_parallel = "none",cluster = NULL, output_frame=FALSE){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  assert_that(all(c("region","year","cases") %in% names(template)))
  #TODO - Set flag if "deaths" in template header
  input_data=input_data_process(input_data,NULL,template)
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0-1")
  assert_that(p_severe_inf >=0.0 && p_severe_inf <=1.0,msg="Severe infection rate must be between 0-1")
  assert_that(p_death_severe_inf >=0.0 && p_death_severe_inf <=1.0,
              msg="Fatality rate of serious infections must be between 0-1")
  assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,msg="Severe infection reporting probability must be between 0-1")
  assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,msg="Fatal infection reporting probability must be between 0-1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  #Prune input data based on regions
  regions=regions_breakdown(template$region)
  input_data=input_data_truncate(input_data,regions)

  #Cross-reference templates with input regions
  xref=template_region_xref(template,input_data$region_labels)

  inv_reps=1/n_reps
  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data
  model_case_values=model_death_values=rep(0,nrow(template))

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(xref$year_data_begin):max(xref$year_end))
    model_output_all=Model_Run_Multi_Input(FOI_spillover = FOI_values,R0 = R0_values,
                                           vacc_data = input_data$vacc_data, pop_data = input_data$pop_data,
                                           years_data = years_data_all, start_SEIRV=start_SEIRV,output_type = "case",
                                           year0 = input_data$years_labels[1],mode_start = mode_start,
                                           vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                           n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(xref$year_data_begin[n_region]:xref$year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV,
                                MoreArgs=list(output_type="case", year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(xref$year_data_begin[n_region]:xref$year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = "case",
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
      #cat("\n\t\tFinished modelling region ",n_region)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile case data
    case_line_list_region=xref$line_list[[n_region]]
    years_case=template$year[case_line_list_region]
    n_lines=length(case_line_list_region)

    for(n_rep in 1:n_reps){
      rep_cases=rep_deaths=rep(0,n_lines)
      for(n_line in 1:n_lines){
        pts=c(1:t_pts)[model_output$year==years_case[n_line]]
        infs=sum(model_output$C[n_rep,pts])
        if(deterministic){
          severe_infs=floor(infs)*p_severe_inf
          deaths=severe_infs*p_death_severe_inf
          rep_cases[n_line]=rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)
        } else {
          severe_infs=rbinom(1,floor(infs),p_severe_inf)
          deaths=rbinom(1,severe_infs,p_death_severe_inf)
          rep_deaths[n_line]=rbinom(1,deaths,p_rep_death)
          rep_cases[n_line]=rep_deaths[n_line]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
        }
      }

      model_case_values[case_line_list_region]=model_case_values[case_line_list_region]+rep_cases
      model_death_values[case_line_list_region]=model_death_values[case_line_list_region]+rep_deaths
    }
  }

  model_case_values=model_case_values*inv_reps
  model_death_values=model_death_values*inv_reps

  if(output_frame) { #Output complete frame of data
    return(data.frame(region=template$region,year=template$year,
                      cases=model_case_values,deaths=model_death_values))
  } else { #Minimal output for MCMC
    return(list(model_case_values=model_case_values,model_death_values=model_death_values))
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
#' @param input_data List of population and vaccination data for multiple regions
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param template Burden data in VIMC format, with regions, years, minimum and maximum age, and life expectancy for each line
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param YLD_per_case TBA
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' '
#' @export
#'
Generate_VIMC_Burden_Dataset <- function(input_data = list(), FOI_values = c(), R0_values = c(), template = NULL,
                                         vaccine_efficacy = 1.0, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                         YLD_per_case = 0.006486, mode_start = 1, start_SEIRV = NULL,
                                         dt = 1.0, n_reps = 1, deterministic = FALSE,mode_parallel = "none",cluster = NULL){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  assert_that(all(c("region","year","age_min","age_max","life_exp") %in% names(template)))
  input_data=input_data_process(input_data,NULL,template)
  assert_that(vaccine_efficacy>=0.0 && vaccine_efficacy<=1.0,msg="Vaccine efficacy must be between 0-1")
  assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
  assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,msg="Fatality rate of severe infections must be between 0-1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  #Prune input data based on regions
  regions=regions_breakdown(template$region)
  input_data=input_data_truncate(input_data,regions)

  #Cross-reference templates with input regions
  xref=template_region_xref(template,input_data$region_labels)

  inv_reps=1/n_reps
  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data
  model_case_values=model_death_values=model_dalys_values=rep(0,nrow(template))

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(xref$year_data_begin):max(xref$year_end))
    model_output_all=Model_Run_Multi_Input(FOI_spillover = FOI_values,R0 = R0_values,
                                           vacc_data = input_data$vacc_data, pop_data = input_data$pop_data,
                                           years_data = years_data_all, start_SEIRV=start_SEIRV,output_type = "case_alt",
                                           year0 = input_data$years_labels[1],mode_start = mode_start,
                                           vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                           n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(xref$year_data_begin[n_region]:xref$year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV,
                                MoreArgs=list(output_type="case_alt", year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(xref$year_data_begin[n_region]:xref$year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = "case_alt",
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
      #cat("\n\t\tFinished modelling region ",n_region)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts_all=c(1:length(model_output$year))

    #Compile case data
    case_line_list_region=xref$line_list[[n_region]]
    years_case=template$year[case_line_list_region]
    n_lines=length(case_line_list_region)

    for(n_rep in 1:n_reps){
      cases=deaths=rep(0,n_lines)
      for(n_line in 1:n_lines){
        line=case_line_list_region[n_line]
        year=years_case[n_line]
        t_pts=t_pts_all[model_output$year==years_case[n_line]]
        n_age_min=findInterval(template$age_min[line],input_data$age_labels)
        n_age_max=findInterval(template$age_max[line],input_data$age_labels)
        age_pts=c(n_age_min:n_age_max)
        infs=sum(model_output$C[age_pts,n_rep,t_pts])
        if(deterministic){
          cases[n_line]=floor(infs)*p_severe_inf
          deaths[n_line]=cases[n_line]*p_death_severe_inf
        } else {
          cases[n_line]=rbinom(1,floor(infs),p_severe_inf)
          deaths[n_line]=rbinom(1,cases[n_line],p_death_severe_inf)
        }
      }

      model_case_values[case_line_list_region]=model_case_values[case_line_list_region]+cases
      model_death_values[case_line_list_region]=model_death_values[case_line_list_region]+deaths
    }
  }

  model_case_values=model_case_values*inv_reps
  model_death_values=model_death_values*inv_reps
  model_YLL_values=model_death_values*template$life_exp
  model_dalys_values=(model_case_values*YLD_per_case)+model_YLL_values

  return(data.frame(region=template$region,year=template$year,age_min=template$age_min,age_max=template$age_max,
                    cases=model_case_values,dalys=model_dalys_values,deaths=model_death_values,YLL=model_YLL_values))
}
#-------------------------------------------------------------------------------
#' @title Generate_Multiple Datasets
#'
#' @description Generate multiple datasets for multiple sets of input parameter values
#'
#' @details This function is used to generate annual serological and/or case/death data based on templates for
#'   multiple sets of parameters (spillover FOI, R0, case/death reporting probabilities, vaccine efficacy)
#'
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param FOI_values Array of sets of values of the force of infection due to spillover from sylvatic reservoir;
#'   dimensions [TBA]
#' @param R0_values Array of values of the basic reproduction number for human-human transmission;
#'   dimensions [TBA]
#' @param sero_template Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#'   and number of samples
#' @param case_template Annual reported case/death data template - data frame with region and year
#' @param vaccine_efficacy Vector of values of vaccine efficacy [TBA]
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Vector of values of probability of reporting of a severe but non-fatal infection [TBA]
#' @param p_rep_death Vetor of value of probability of reporting of a fatal infection [TBA]
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel Set mode for parallelization, if any:
#'   If mode_parallel="none", no parallelization
#'   If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
#'   If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types
#' @param cluster Cluster of threads to use if mode_parallel="clusterMap"
#' '
#' @export
#'
Generate_Multiple_Datasets <- function(input_data = list(),FOI_values = list(),R0_values = list(),sero_template = NULL,case_template = NULL,
                                       vaccine_efficacy = c(1.0), p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = c(1.0),
                                       p_rep_death = c(1.0),mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                                       mode_parallel = "none",cluster = NULL){

  assert_that(length(dim(FOI_values))==2)
  assert_that(length(dim(FOI_values))==2)
  n_param_sets=dim(FOI_values)[2]
  assert_that(dim(R0_values)[2]==n_param_sets)
  assert_that(length(vaccine_efficacy)==n_param_sets)
  assert_that(length(p_rep_severe)==n_param_sets)
  assert_that(length(p_rep_death)==n_param_sets)

  if(is.null(sero_template)==FALSE){model_sero_values=array(NA,dim=c(nrow(sero_template),n_param_sets))}else(model_sero_values=NULL)
  if(is.null(case_template)==FALSE){model_case_values=model_death_values=array(NA,dim=c(nrow(case_template),n_param_sets))}else{
    model_case_values=model_death_values=NULL
  }
  cat("\nRunning set:\n")
  for(set in 1:n_param_sets){
    cat("\t",set,sep="")
    dataset_single <- Generate_Dataset(input_data,FOI_values[,set],R0_values[,set],sero_template,case_template,vaccine_efficacy[set],
                                       p_severe_inf,p_death_severe_inf,p_rep_severe[set],p_rep_death[set],mode_start,start_SEIRV,dt,n_reps,
                                       deterministic,mode_parallel,cluster,output_frame=FALSE)
    if(is.null(sero_template)==FALSE){
      model_sero_values[,set]=dataset_single$model_sero_values
    }
    if(is.null(case_template)==FALSE){
      model_case_values[,set]=dataset_single$model_case_values
      model_death_values[,set]=dataset_single$model_death_values
    }
  }
  cat("\nDatasets complete.\n")

  return(list(model_sero_values=model_sero_values,model_case_values=model_case_values,
              model_death_values=model_death_values))
}
