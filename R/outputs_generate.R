# R file for functions used to generate sets of annual serological data, annual case/death data
# and annual burden (VIMC format) data
#-------------------------------------------------------------------------------
#' @title Generate_Dataset
#'
#' @description Generate serological and/or annual case/death data
#'
#' @details This function is used to generate serological and/or annual case/death based on observed or dummy data sets;
#' it is normally used by the single_like_calc() and data_match_single() functions. The separate Generate_Sero_dataset
#' and Generate_Case_Dataset functions can be used when only one or the other type is required; this function exists
#' to cover instances where serological and case data may need to be generated for the same region.
#'
#' @param input_data List of population and vaccination data for multiple regions
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
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
#' '
#' @export
#'
Generate_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),obs_sero_data = NULL,obs_case_data = NULL,
                             vaccine_efficacy = 1.0, p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = 1.0,p_rep_death = 1.0,
                             mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE, mode_parallel = "none",
                             cluster = NULL){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  if(is.null(input_data$flag_sero)){input_data=input_data_process(input_data,obs_sero_data,obs_case_data)}
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE),
              msg="Need at least one of obs_sero_data or obs_case_data")
  if(is.null(obs_sero_data)==FALSE){
    assert_that(all(c("age_min","age_max","samples","positives","vc_factor","region") %in% names(obs_sero_data)))
  }
  if(is.null(obs_case_data)==FALSE){
    assert_that(all(c("region","year","cases") %in% names(obs_case_data)))
    #TODO - Flag if "deaths" in header
    assert_that(p_severe_inf >=0.0 && p_severe_inf <=1.0,msg="Severe infection rate must be between 0 and 1")
    assert_that(p_death_severe_inf >=0.0 && p_death_severe_inf <=1.0,msg="Fatality rate of serious infections must be between 0 and 1")
    assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,msg="Severe infection reporting probability must be between 0 and 1")
    assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,msg="Fatal infection reporting probability must be between 0 and 1")
  }
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  inv_reps=1/n_reps
  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)){model_sero_data=NULL} else {
    blank=rep(0,nrow(obs_sero_data))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  }
  if(is.null(obs_case_data)){model_case_values=model_death_values=NA} else {
    model_case_values=model_death_values=rep(0,nrow(obs_case_data))
  }

  #Set up vector of output types to get from model if needed
  if(mode_parallel %in% c("none","clusterMap")){
    output_types=rep(NA,n_regions)
    for(n_region in 1:n_regions){
      if(input_data$flag_case[n_region]==1){
        if(input_data$flag_sero[n_region]==1){output_types[n_region]="case+sero"} else{output_types[n_region]="case"}
      } else {output_types[n_region]="sero"}
    }
  }

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(input_data$year_data_begin):max(input_data$year_end))
    if(any(input_data$flag_sero==1)){if(any(input_data$flag_case==1)){output_type="case+sero"} else{
      output_type="sero"}} else {output_type="case"}
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
      years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
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
                               years_data = c(input_data$year_data_begin[n_region]:input_data$year_end[n_region]),
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
    if(input_data$flag_case[n_region]==1){
      case_line_list=input_data$case_line_list[[n_region]]
      years_outbreak=obs_case_data$year[case_line_list]
      n_years_outbreak=length(case_line_list)

      for(n_rep in 1:n_reps){
        rep_cases=rep_deaths=rep(0,n_years_outbreak)
        for(n_year in 1:n_years_outbreak){
          pts=c(1:t_pts)[model_output$year==years_outbreak[n_year]]
          infs=sum(model_output$C[n_rep,pts])
          if(deterministic){
            severe_infs=floor(infs)*p_severe_inf
            deaths=severe_infs*p_death_severe_inf
            rep_deaths[n_year]=round(deaths*p_rep_death)
            rep_cases[n_year]=rep_deaths[n_year]+round((severe_infs-deaths)*p_rep_severe)

          } else {
            severe_infs=rbinom(1,floor(infs),p_severe_inf)
            deaths=rbinom(1,severe_infs,p_death_severe_inf)
            rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
            rep_cases[n_year]=rep_deaths[n_year]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
          }
        }

        model_case_values[case_line_list]=model_case_values[case_line_list]+rep_cases
        model_death_values[case_line_list]=model_death_values[case_line_list]+rep_deaths
      }
    }

    #Compile seroprevalence data if necessary
    if(input_data$flag_sero[n_region]==1){
      sero_line_list=input_data$sero_line_list[[n_region]]
      for(n_rep in 1:n_reps){
        sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_output,n_rep)
        model_sero_data$samples[sero_line_list]=model_sero_data$samples[sero_line_list]+sero_results$samples
        model_sero_data$positives[sero_line_list]=model_sero_data$positives[sero_line_list]+sero_results$positives
      }
    }
    # model_output<-NULL
    # gc()
  }
  # if(mode_parallel!="none"){
  #   model_output_all<-NULL
  #   gc()
  # }

  if(any(input_data$flag_sero>0)){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(any(input_data$flag_case>0)){
    model_case_values=model_case_values*inv_reps
    model_death_values=model_death_values*inv_reps
  }

  # if(output_flag==0){ #Minimal output for MCMC
    return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
                model_death_values=model_death_values))
  # } else { #TBA - Option for outputting complete frame of data
  #   return(list(model_sero_data=data.frame(region=obs_sero_data$region,year=obs_sero_data$year,
  #                                          age_min=obs_sero_data$age_min,age_max=obs_sero_data$age_max,
  #                                          samples=obs_sero_data$samples,positives=obs_sero_data$samples*model_sero_data$sero),
  #               model_case_data=data.frame(region=obs_case_data$region,year=obs_case_data$year,
  #                                          cases=model_case_values,deaths=model_death_values)))
  # }
}
#-------------------------------------------------------------------------------
#' @title Generate_Sero_Dataset
#'
#' @description Generate serological dataset
#'
#' @details This function is used to generate serological, data for multiple regions/parameter sets based
#'   on a template.
#'
#' @param input_data List of population and vaccination data for multiple regions
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
#' '
#' @export
#'
Generate_Sero_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),template = NULL,
                                  vaccine_efficacy = 1.0, mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1,
                                  deterministic = FALSE, mode_parallel = "none",cluster = NULL){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  assert_that(all(c("age_min","age_max","samples","positives","vc_factor","region") %in% names(template)))
  input_data=input_data_process(input_data,template,NULL)
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

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
    years_data_all=c(min(input_data$year_data_begin):max(input_data$year_end))
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
      years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
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
                               years_data = c(input_data$year_data_begin[n_region]:input_data$year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = "sero",
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile seroprevalence data
    sero_line_list=input_data$sero_line_list[[n_region]]
    for(n_rep in 1:n_reps){
      sero_results=sero_calculate2(template[sero_line_list,],model_output,n_rep)
      model_sero_data$samples[sero_line_list]=model_sero_data$samples[sero_line_list]+sero_results$samples
      model_sero_data$positives[sero_line_list]=model_sero_data$positives[sero_line_list]+sero_results$positives
    }
    # model_output<-NULL
    # gc()
  }
  # if(mode_parallel!="none"){
  #   model_output_all<-NULL
  #   gc()
  # }
  model_sero_data$sero=model_sero_data$positives/model_sero_data$samples

  # if(output_flag==0){ #Minimal output for MCMC
    return(model_sero_data)
  # } else { #TBA - Option for outputting complete frame of data
  #   return(data.frame(region=template$region,year=template$year,
  #                     age_min=template$age_min,age_max=template$age_max,
  #                     samples=template$samples,positives=template$samples*model_sero_data$sero))
  # }
}
#-------------------------------------------------------------------------------
#' @title Generate_Case_Dataset
#'
#' @description Generate annual case/death data
#'
#' @details This function is used to generate annual case/death data based on a template
#'
#' @param input_data List of population and vaccination data for multiple regions
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
#' '
#' @export
#'
Generate_Case_Dataset <- function(input_data = list(),FOI_values = c(),R0_values = c(),template = NULL,vaccine_efficacy = 1.0,
                                  p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = 1.0,p_rep_death = 1.0,
                                  mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                                  mode_parallel = "none",cluster = NULL){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see [TBA] )"))
  assert_that(all(c("region","year","cases") %in% names(template)))
  #TODO - Set flag if "deaths" in template header
  input_data=input_data_process(input_data,NULL,template)
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  assert_that(p_severe_inf >=0.0 && p_severe_inf <=1.0,msg="Severe infection rate must be between 0 and 1")
  assert_that(p_death_severe_inf >=0.0 && p_death_severe_inf <=1.0,
              msg="Fatality rate of serious infections must be between 0 and 1")
  assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,msg="Severe infection reporting probability must be between 0 and 1")
  assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,msg="Fatal infection reporting probability must be between 0 and 1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

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
    years_data_all=c(min(input_data$year_data_begin):max(input_data$year_end))
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
      years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
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
                               years_data = c(input_data$year_data_begin[n_region]:input_data$year_end[n_region]),
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
    case_line_list=input_data$case_line_list[[n_region]]
    years_case=template$year[case_line_list]
    n_years_case=length(case_line_list)

    for(n_rep in 1:n_reps){
      rep_cases=rep_deaths=rep(0,n_years_case)
      for(n_year in 1:n_years_case){
        pts=c(1:t_pts)[model_output$year==years_case[n_year]]
        infs=sum(model_output$C[n_rep,pts])
        if(deterministic){
          severe_infs=floor(infs)*p_severe_inf
          deaths=severe_infs*p_death_severe_inf
          rep_deaths[n_year]=round(deaths*p_rep_death)
          rep_cases[n_year]=rep_deaths[n_year]+round((severe_infs-deaths)*p_rep_severe)

        } else {
          severe_infs=rbinom(1,floor(infs),p_severe_inf)
          deaths=rbinom(1,severe_infs,p_death_severe_inf)
          rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
          rep_cases[n_year]=rep_deaths[n_year]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
        }
      }

      model_case_values[case_line_list]=model_case_values[case_line_list]+rep_cases
      model_death_values[case_line_list]=model_death_values[case_line_list]+rep_deaths
    }

    # model_output<-NULL
    # gc()
  }
  # if(mode_parallel!="none"){
  #   model_output_all<-NULL
  #   gc()
  # }
  model_case_values=model_case_values*inv_reps
  model_death_values=model_death_values*inv_reps

  # if(output_flag==0){ #Minimal output for MCMC
    return(list(model_case_values=model_case_values,model_death_values=model_death_values))
  # } else { #TBA - Option for outputting complete frame of data
  #   return(data.frame(region=template$region,year=template$year,
  #                     cases=model_case_values,deaths=model_death_values))
  # }

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
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  assert_that(p_severe_inf >=0.0 && p_severe_inf <=1.0,msg="Severe infection rate must be between 0 and 1")
  assert_that(p_death_severe_inf >=0.0 && p_death_severe_inf <=1.0,msg="Fatality rate of serious infections must be between 0 and 1")
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

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
    years_data_all=c(min(input_data$year_data_begin):max(input_data$year_end))
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
      years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
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
                               years_data = c(input_data$year_data_begin[n_region]:input_data$year_end[n_region]),
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
    case_line_list=input_data$case_line_list[[n_region]]
    years_case=template$year[case_line_list]
    n_lines=length(case_line_list)

    for(n_rep in 1:n_reps){
      cases=dalys=deaths=rep(0,n_lines)
      for(n_line in 1:n_lines){
        line=case_line_list[n_line]
        year=years_case[n_line]
        t_pts=t_pts_all[model_output$year==years_case[n_line]]
        n_age_min=match(template$age_min[line],input_data$age_labels)
        n_age_max=match(template$age_max[line],input_data$age_labels)
        age_pts=c(n_age_min:n_age_max) #TBC - Currently just uses all ages
        infs=sum(model_output$C[age_pts,n_rep,t_pts])
        if(deterministic){
          cases[n_line]=floor(infs)*p_severe_inf
          deaths[n_line]=cases[n_line]*p_death_severe_inf
        } else {
          cases[n_line]=rbinom(1,floor(infs),p_severe_inf)
          deaths[n_line]=rbinom(1,cases[n_line],p_death_severe_inf)
        }
        dalys[n_line]=(cases[n_line]*YLD_per_case)
      }

      model_case_values[case_line_list]=model_case_values[case_line_list]+cases
      model_death_values[case_line_list]=model_death_values[case_line_list]+deaths
    }

    # model_output<-NULL
    # gc()
  }
  # if(mode_parallel!="none"){
  #   model_output_all<-NULL
  #   gc()
  # }
  model_case_values=model_case_values*inv_reps
  model_death_values=model_death_values*inv_reps
  model_dalys_values=(model_case_values*YLD_per_case)+(model_death_values*template$life_exp)

  # if(output_flag==0){ #Minimal output for MCMC
    return(list(cases=model_case_values,dalys=model_dalys_values,deaths=model_death_values))
  # } else { #TBA - Option for outputting complete frame of data
  #   return(data.frame(region=template$region,year=template$year,age_min=template$age_min,age_max=template$age_max,
  #                     cases=model_case_values,dalys=model_dalys_values,deaths=model_death_values))
  # }
}
